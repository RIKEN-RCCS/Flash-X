#!/usr/bin/env python3
"""Expand ``@name(...)@`` macros using section-based definition files.

Definition files have this form::

    [spc]
    args=ii,jj,kk
    definition=
      ii,jj,kk

Definitions may invoke other macros.  Formal arguments can be written as
standalone identifiers, as ``${name}``, or concatenated with other formal
arguments.  For example, with formals ``ib,io``, the token ``ibio`` expands
to the concatenation of the two actual arguments.
"""

from __future__ import annotations

import argparse
import os
import re
import stat
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


IDENTIFIER = r"[A-Za-z_][A-Za-z0-9_]*"
IDENTIFIER_RE = re.compile(rf"^{IDENTIFIER}$")
SECTION_RE = re.compile(rf"^\s*\[({IDENTIFIER})\]\s*$")
FIELD_RE = re.compile(r"^\s*(args|definition)\s*=\s*(.*)$")
SUBSTITUTION_RE = re.compile(
    rf"\$\{{({IDENTIFIER})\}}|({IDENTIFIER})"
)


class MacroError(Exception):
    """Base class for definition and expansion errors."""


class DefinitionError(MacroError):
    """A macro-definition file is malformed."""


class ExpansionError(MacroError):
    """A source macro cannot be expanded safely."""


@dataclass(frozen=True)
class Macro:
    name: str
    arguments: tuple[str, ...]
    definition: str
    path: Path
    line: int


def _strip_common_indent(lines: list[str]) -> str:
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()

    nonblank = [line for line in lines if line.strip()]
    if not nonblank:
        return ""

    indent = min(len(line) - len(line.lstrip()) for line in nonblank)
    return "\n".join(line[indent:] if line.strip() else "" for line in lines)


def _parse_section(
    path: Path,
    name: str,
    section_line: int,
    body: list[tuple[int, str]],
) -> Macro:
    arguments: tuple[str, ...] = ()
    saw_arguments = False
    definition_lines: list[str] | None = None

    for line_number, line in body:
        if definition_lines is not None:
            definition_lines.append(line)
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue

        field = FIELD_RE.match(line)
        if field is None:
            raise DefinitionError(
                f"{path}:{line_number}: expected args= or definition= in [{name}]"
            )

        key, value = field.groups()
        if key == "args":
            if saw_arguments:
                raise DefinitionError(
                    f"{path}:{line_number}: duplicate args= in [{name}]"
                )
            saw_arguments = True
            parsed = tuple(part.strip() for part in value.split(",") if part.strip())
            if value.strip() and len(parsed) != len(value.split(",")):
                raise DefinitionError(
                    f"{path}:{line_number}: empty formal argument in [{name}]"
                )
            invalid = [argument for argument in parsed if not IDENTIFIER_RE.fullmatch(argument)]
            if invalid:
                raise DefinitionError(
                    f"{path}:{line_number}: invalid formal argument {invalid[0]!r}"
                )
            if len(set(parsed)) != len(parsed):
                raise DefinitionError(
                    f"{path}:{line_number}: duplicate formal argument in [{name}]"
                )
            arguments = parsed
        else:
            definition_lines = [value] if value else []

    if definition_lines is None:
        raise DefinitionError(
            f"{path}:{section_line}: macro [{name}] has no definition="
        )

    return Macro(
        name=name,
        arguments=arguments,
        definition=_strip_common_indent(definition_lines),
        path=path,
        line=section_line,
    )


def parse_definition_file(path: Path) -> list[Macro]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DefinitionError(f"cannot read {path}: {exc}") from exc

    macros: list[Macro] = []
    current_name: str | None = None
    current_line = 0
    current_body: list[tuple[int, str]] = []

    for line_number, line in enumerate(lines, start=1):
        section = SECTION_RE.match(line)
        if section:
            if current_name is not None:
                macros.append(
                    _parse_section(
                        path, current_name, current_line, current_body
                    )
                )
            current_name = section.group(1)
            current_line = line_number
            current_body = []
        elif current_name is None:
            stripped = line.strip()
            if stripped and not stripped.startswith(("#", ";")):
                raise DefinitionError(
                    f"{path}:{line_number}: content outside a [macro] section"
                )
        else:
            current_body.append((line_number, line))

    if current_name is not None:
        macros.append(_parse_section(path, current_name, current_line, current_body))

    if not macros:
        raise DefinitionError(f"{path}: no macro definitions found")
    return macros


def load_definitions(paths: Iterable[Path]) -> dict[str, Macro]:
    macros: dict[str, Macro] = {}
    for path in paths:
        for macro in parse_definition_file(path):
            if macro.name in macros:
                previous = macros[macro.name]
                raise DefinitionError(
                    f"duplicate macro [{macro.name}] at {macro.path}:{macro.line}; "
                    f"previously defined at {previous.path}:{previous.line}"
                )
            macros[macro.name] = macro
    return macros


def split_arguments(text: str) -> list[str]:
    if not text.strip():
        return []

    arguments: list[str] = []
    start = 0
    delimiters: list[str] = []
    quote: str | None = None
    macro_open = False
    matching = {")": "(", "]": "[", "}": "{"}

    i = 0
    while i < len(text):
        char = text[i]

        if quote is not None:
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2
                    continue
                quote = None
            i += 1
            continue

        if char in ("'", '"'):
            quote = char
        elif char == "@":
            macro_open = not macro_open
        elif not macro_open:
            if char in "([{":
                delimiters.append(char)
            elif char in ")]}":
                if not delimiters or delimiters.pop() != matching[char]:
                    raise ExpansionError(f"unbalanced delimiter in argument list {text!r}")
            elif char == "," and not delimiters:
                arguments.append(text[start:i].strip())
                start = i + 1
        i += 1

    if quote is not None or macro_open or delimiters:
        raise ExpansionError(f"unbalanced argument list {text!r}")

    arguments.append(text[start:].strip())
    if any(not argument for argument in arguments):
        raise ExpansionError(f"empty argument in {text!r}")
    return arguments


def parse_invocation(text: str) -> tuple[str, list[str] | None]:
    text = text.strip()
    if IDENTIFIER_RE.fullmatch(text):
        return text, None

    open_parenthesis = text.find("(")
    if open_parenthesis < 1 or not text.endswith(")"):
        raise ExpansionError(f"malformed macro invocation @{text}@")

    name = text[:open_parenthesis].strip()
    if not IDENTIFIER_RE.fullmatch(name):
        raise ExpansionError(f"invalid macro name in @{text}@")
    return name, split_arguments(text[open_parenthesis + 1 : -1])


def find_closing_at(text: str, opening: int) -> int:
    delimiters: list[str] = []
    quote: str | None = None
    matching = {")": "(", "]": "[", "}": "{"}
    i = opening + 1

    while i < len(text):
        char = text[i]
        if quote is not None:
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2
                    continue
                quote = None
        elif char in ("'", '"'):
            quote = char
        elif char in "([{":
            delimiters.append(char)
        elif char in ")]}":
            if not delimiters or delimiters.pop() != matching[char]:
                raise ExpansionError(
                    f"unbalanced delimiter near character {i}"
                )
        elif char == "@":
            if not delimiters:
                return i
            i = find_closing_at(text, i)
        i += 1

    raise ExpansionError(f"opening '@' at character {opening} is not closed")


def _segment_formals(token: str, names: tuple[str, ...]) -> list[str] | None:
    memo: dict[int, list[list[str]]] = {}

    def visit(position: int) -> list[list[str]]:
        if position == len(token):
            return [[]]
        if position in memo:
            return memo[position]

        results: list[list[str]] = []
        for name in names:
            if token.startswith(name, position):
                results.extend([[name, *tail] for tail in visit(position + len(name))])
                if len(results) > 1:
                    break
        memo[position] = results
        return results

    segmentations = visit(0)
    if not segmentations:
        return None
    if len(segmentations) > 1:
        raise ExpansionError(
            f"ambiguous concatenated formal arguments in token {token!r}; "
            "use ${argument} syntax"
        )
    return segmentations[0]


def substitute_arguments(macro: Macro, actual: list[str]) -> str:
    values = dict(zip(macro.arguments, actual))

    def replace(match: re.Match[str]) -> str:
        explicit_name = match.group(1)
        token = match.group(2)
        if explicit_name is not None:
            if explicit_name not in values:
                raise ExpansionError(
                    f"macro [{macro.name}] references unknown formal "
                    f"argument {explicit_name!r}"
                )
            return values[explicit_name]

        assert token is not None
        if token in values:
            return values[token]
        segments = _segment_formals(token, macro.arguments)
        if segments is None or len(segments) < 2:
            return token
        return "".join(values[name] for name in segments)

    return SUBSTITUTION_RE.sub(replace, macro.definition)


class MacroExpander:
    def __init__(self, macros: dict[str, Macro], max_depth: int = 100):
        self.macros = macros
        self.max_depth = max_depth
        self.expansion_count = 0

    def expand_text(self, text: str, stack: tuple[str, ...] = ()) -> str:
        output: list[str] = []
        cursor = 0
        while cursor < len(text):
            opening = text.find("@", cursor)
            if opening == -1:
                output.append(text[cursor:])
                break
            output.append(text[cursor:opening])
            closing = find_closing_at(text, opening)
            output.append(self.expand_invocation(text[opening + 1 : closing], stack))
            cursor = closing + 1
        return "".join(output)

    def expand_invocation(self, invocation: str, stack: tuple[str, ...]) -> str:
        name, supplied = parse_invocation(invocation)
        if name not in self.macros:
            raise ExpansionError(f"undefined macro {name!r}")
        if name in stack:
            raise ExpansionError(
                "recursive macro expansion: " + " -> ".join((*stack, name))
            )
        if len(stack) >= self.max_depth:
            raise ExpansionError(
                f"maximum expansion depth {self.max_depth} exceeded"
            )

        macro = self.macros[name]
        actual = [] if supplied is None else supplied
        if len(actual) != len(macro.arguments):
            raise ExpansionError(
                f"macro {name!r} expects {len(macro.arguments)} argument(s), "
                f"received {len(actual)}"
            )

        new_stack = (*stack, name)
        expanded_actual = [self.expand_text(item, new_stack) for item in actual]
        substituted = substitute_arguments(macro, expanded_actual)
        self.expansion_count += 1
        return self.expand_text(substituted, new_stack)


def write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        output_mode = stat.S_IMODE(path.stat().st_mode)
    except FileNotFoundError:
        output_mode = 0o644
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary.write(text)
            temporary_name = temporary.name
        os.chmod(temporary_name, output_mode)
        os.replace(temporary_name, path)
    except OSError as exc:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass
        raise ExpansionError(f"cannot write {path}: {exc}") from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Expand @name(...)@ macros in a text source file."
    )
    parser.add_argument(
        "-d", "--definitions", nargs="+", required=True, metavar="FILE"
    )
    parser.add_argument("-i", "--input", required=True, metavar="FILE")
    parser.add_argument("-o", "--output", required=True, metavar="FILE")
    parser.add_argument("--max-depth", type=int, default=100, metavar="N")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    if options.max_depth < 1:
        parser.error("--max-depth must be at least 1")

    input_path = Path(options.input)
    output_path = Path(options.output)
    try:
        macros = load_definitions(Path(item) for item in options.definitions)
        try:
            source = input_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise ExpansionError(f"cannot read {input_path}: {exc}") from exc

        expander = MacroExpander(macros, max_depth=options.max_depth)
        expanded = expander.expand_text(source)
        if "@" in expanded:
            raise ExpansionError("unexpanded '@' delimiter remains in output")
        write_atomic(output_path, expanded)
    except MacroError as exc:
        print(f"{Path(sys.argv[0]).name}: error: {exc}", file=sys.stderr)
        return 1

    print(
        f"expanded {expander.expansion_count} invocation(s) from "
        f"{input_path} into {output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
