#!/usr/bin/env python3
"""
Embedded macro expander.

Macro invocation syntax:
    @macro_name@
    @macro_name(arg1,arg2)@

Example definition file:

    space_index
    definition=
       i,j,k

    index5d
    args=var,blk
    definition=
       var,@space_index@,blk

Example source:

    array1(@index5d(rho,me)@,@space_index@) = array2(@space_index@)

Usage:

    python macro_expand.py \
        --definitions general.mac project.mac \
        --input input.F90 \
        --output output.F90

Later definition files override earlier ones by default.
Use --reject-duplicates to make duplicate macro names an error.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
IDENTIFIER_TOKEN_RE = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")


class MacroError(Exception):
    """Base class for all macro-expansion errors."""


class DefinitionError(MacroError):
    """Raised for malformed macro-definition files."""


class ExpansionError(MacroError):
    """Raised for malformed macro invocations or expansion failures."""


@dataclass(frozen=True)
class Macro:
    name: str
    args: tuple[str, ...]
    definition: str
    source_file: Path
    source_line: int


def strip_common_indent(lines: list[str]) -> str:
    """
    Remove surrounding blank lines and common leading indentation.

    Internal line breaks are preserved.
    """
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()

    if not lines:
        return ""

    nonblank = [line for line in lines if line.strip()]
    indent = min(len(line) - len(line.lstrip()) for line in nonblank)

    return "\n".join(
        line[indent:] if line.strip() else ""
        for line in lines
    )


def parse_macro_file(path: Path) -> list[Macro]:
    """
    Parse one macro-definition file.

    Supported structure:

        macro_name
        args=a,b
        definition=
           replacement text

    The args= line is optional. A definition may start on the same line:

        definition=i,j,k

    A definition continues until the next non-indented line that is not blank
    and is not part of the current macro's fields.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise DefinitionError(f"Cannot read definition file {path}: {exc}") from exc

    lines = text.splitlines()
    macros: list[Macro] = []
    i = 0

    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()

        if not stripped or stripped.startswith("#"):
            i += 1
            continue

        if raw != raw.lstrip():
            raise DefinitionError(
                f"{path}:{i + 1}: expected a macro name at the left margin"
            )

        name = stripped
        if not IDENTIFIER_RE.fullmatch(name):
            raise DefinitionError(
                f"{path}:{i + 1}: invalid macro name {name!r}"
            )

        macro_line = i + 1
        i += 1
        args: tuple[str, ...] = ()
        definition_lines: list[str] | None = None

        while i < len(lines):
            raw = lines[i]
            stripped = raw.strip()

            if not stripped or stripped.startswith("#"):
                i += 1
                continue

            if raw == raw.lstrip() and not (
                stripped.startswith("args=")
                or stripped.startswith("definition=")
            ):
                break

            if stripped.startswith("args="):
                if definition_lines is not None:
                    raise DefinitionError(
                        f"{path}:{i + 1}: args= must appear before definition="
                    )
                arg_text = stripped[len("args="):].strip()
                if not arg_text:
                    args = ()
                else:
                    parsed_args = [part.strip() for part in arg_text.split(",")]
                    for arg in parsed_args:
                        if not IDENTIFIER_RE.fullmatch(arg):
                            raise DefinitionError(
                                f"{path}:{i + 1}: invalid formal argument {arg!r}"
                            )
                    if len(set(parsed_args)) != len(parsed_args):
                        raise DefinitionError(
                            f"{path}:{i + 1}: duplicate formal argument"
                        )
                    args = tuple(parsed_args)
                i += 1
                continue

            if stripped.startswith("definition="):
                if definition_lines is not None:
                    raise DefinitionError(
                        f"{path}:{i + 1}: duplicate definition= field"
                    )

                after_equals = raw[raw.find("definition=") + len("definition="):]
                definition_lines = []
                if after_equals.strip():
                    definition_lines.append(after_equals.lstrip())

                i += 1
                while i < len(lines):
                    continuation = lines[i]

                    if not continuation.strip():
                        definition_lines.append("")
                        i += 1
                        continue

                    if continuation == continuation.lstrip():
                        break

                    definition_lines.append(continuation)
                    i += 1
                continue

            raise DefinitionError(
                f"{path}:{i + 1}: expected args= or definition="
            )

        if definition_lines is None:
            raise DefinitionError(
                f"{path}:{macro_line}: macro {name!r} has no definition="
            )

        macros.append(
            Macro(
                name=name,
                args=args,
                definition=strip_common_indent(definition_lines),
                source_file=path,
                source_line=macro_line,
            )
        )

    return macros


def load_macros(
    definition_paths: Iterable[Path],
    reject_duplicates: bool = False,
) -> dict[str, Macro]:
    """Load and combine macro definitions from multiple files."""
    table: dict[str, Macro] = {}

    for path in definition_paths:
        for macro in parse_macro_file(path):
            if reject_duplicates and macro.name in table:
                previous = table[macro.name]
                raise DefinitionError(
                    f"Duplicate macro {macro.name!r}: "
                    f"{previous.source_file}:{previous.source_line} and "
                    f"{macro.source_file}:{macro.source_line}"
                )
            table[macro.name] = macro

    return table


def split_arguments(text: str) -> list[str]:
    """
    Split a macro argument list on top-level commas.

    Parentheses and nested @...@ invocations are respected.
    """
    if not text.strip():
        return []

    args: list[str] = []
    start = 0
    paren_depth = 0
    macro_depth = 0
    i = 0

    while i < len(text):
        char = text[i]

        if char == "@":
            macro_depth = 1 - macro_depth
            i += 1
            continue

        if macro_depth == 0:
            if char == "(":
                paren_depth += 1
            elif char == ")":
                paren_depth -= 1
                if paren_depth < 0:
                    raise ExpansionError(
                        f"Unbalanced ')' in argument list: {text!r}"
                    )
            elif char == "," and paren_depth == 0:
                args.append(text[start:i].strip())
                start = i + 1

        i += 1

    if macro_depth != 0:
        raise ExpansionError(
            f"Unbalanced '@' in argument list: {text!r}"
        )
    if paren_depth != 0:
        raise ExpansionError(
            f"Unbalanced parentheses in argument list: {text!r}"
        )

    args.append(text[start:].strip())

    if any(arg == "" for arg in args):
        raise ExpansionError(f"Empty macro argument in: {text!r}")

    return args


def parse_invocation(content: str) -> tuple[str, list[str] | None]:
    """
    Parse the text between the outer @ characters.

    Returns:
        (macro_name, None) for @name@
        (macro_name, arguments) for @name(...)@
    """
    content = content.strip()
    if not content:
        raise ExpansionError("Empty macro invocation @@")

    open_paren = content.find("(")
    if open_paren == -1:
        if not IDENTIFIER_RE.fullmatch(content):
            raise ExpansionError(f"Invalid macro invocation @{content}@")
        return content, None

    name = content[:open_paren].strip()
    if not IDENTIFIER_RE.fullmatch(name):
        raise ExpansionError(f"Invalid macro name in @{content}@")

    if not content.endswith(")"):
        raise ExpansionError(
            f"Malformed invocation @{content}@; expected closing ')'"
        )

    arg_text = content[open_paren + 1:-1]
    return name, split_arguments(arg_text)


def find_closing_at(text: str, start: int) -> int:
    """
    Find the closing @ for an invocation beginning at text[start].

    Nested macro invocations inside an argument list are supported.
    """
    paren_depth = 0
    i = start + 1

    while i < len(text):
        char = text[i]

        if char == "(":
            paren_depth += 1
        elif char == ")":
            paren_depth -= 1
            if paren_depth < 0:
                raise ExpansionError(
                    f"Unbalanced ')' near character {i}"
                )
        elif char == "@":
            if paren_depth == 0:
                return i

            nested_end = find_closing_at(text, i)
            i = nested_end

        i += 1

    raise ExpansionError(
        f"Opening '@' at character {start} has no matching closing '@'"
    )


def substitute_arguments(
    definition: str,
    formal_args: tuple[str, ...],
    actual_args: list[str],
) -> str:
    """Replace formal arguments as complete identifiers."""
    mapping = dict(zip(formal_args, actual_args))

    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        return mapping.get(token, token)

    return IDENTIFIER_TOKEN_RE.sub(replace, definition)


class MacroExpander:
    def __init__(self, macros: dict[str, Macro], max_depth: int = 100):
        self.macros = macros
        self.max_depth = max_depth

    def expand_text(
        self,
        text: str,
        expansion_stack: tuple[str, ...] = (),
    ) -> str:
        """Expand every embedded macro invocation in text."""
        output: list[str] = []
        cursor = 0

        while cursor < len(text):
            opening = text.find("@", cursor)
            if opening == -1:
                output.append(text[cursor:])
                break

            output.append(text[cursor:opening])
            closing = find_closing_at(text, opening)
            invocation = text[opening + 1:closing]

            expanded = self.expand_invocation(
                invocation,
                expansion_stack=expansion_stack,
            )
            output.append(expanded)
            cursor = closing + 1

        return "".join(output)

    def expand_invocation(
        self,
        invocation: str,
        expansion_stack: tuple[str, ...],
    ) -> str:
        """Expand one macro invocation."""
        name, supplied_args = parse_invocation(invocation)

        if name not in self.macros:
            raise ExpansionError(f"Undefined macro {name!r}")

        if name in expansion_stack:
            cycle = " -> ".join((*expansion_stack, name))
            raise ExpansionError(f"Recursive macro expansion detected: {cycle}")

        if len(expansion_stack) >= self.max_depth:
            chain = " -> ".join(expansion_stack)
            raise ExpansionError(
                f"Maximum expansion depth {self.max_depth} exceeded: {chain}"
            )

        macro = self.macros[name]

        if macro.args:
            if supplied_args is None:
                raise ExpansionError(
                    f"Macro {name!r} requires {len(macro.args)} argument(s)"
                )
            if len(supplied_args) != len(macro.args):
                raise ExpansionError(
                    f"Macro {name!r} expects {len(macro.args)} argument(s), "
                    f"but received {len(supplied_args)}"
                )
        else:
            if supplied_args is not None and supplied_args:
                raise ExpansionError(
                    f"Macro {name!r} takes no arguments"
                )
            supplied_args = []

        new_stack = (*expansion_stack, name)

        expanded_actual_args = [
            self.expand_text(arg, expansion_stack=new_stack)
            for arg in supplied_args
        ]

        substituted = substitute_arguments(
            macro.definition,
            macro.args,
            expanded_actual_args,
        )

        return self.expand_text(
            substituted,
            expansion_stack=new_stack,
        )


def expand_file(
    input_path: Path,
    output_path: Path,
    expander: MacroExpander,
) -> None:
    """Read, expand, and write one source file."""
    try:
        source = input_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ExpansionError(f"Cannot read input file {input_path}: {exc}") from exc

    expanded = expander.expand_text(source)

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(expanded, encoding="utf-8")
    except OSError as exc:
        raise ExpansionError(
            f"Cannot write output file {output_path}: {exc}"
        ) from exc


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Expand embedded @macro@ and @macro(args)@ invocations."
    )
    parser.add_argument(
        "-d",
        "--definitions",
        nargs="+",
        required=True,
        metavar="FILE",
        help=(
            "Macro-definition files. Later files override earlier files "
            "unless --reject-duplicates is used."
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        help="Input source file containing embedded macros.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="FILE",
        help="Output path for the expanded source.",
    )
    parser.add_argument(
        "--reject-duplicates",
        action="store_true",
        help="Treat duplicate macro definitions as an error.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=100,
        help="Maximum recursive expansion depth; default: 100.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    options = parser.parse_args(argv)

    if options.max_depth < 1:
        parser.error("--max-depth must be at least 1")

    definition_paths = [Path(path) for path in options.definitions]
    input_path = Path(options.input)
    output_path = Path(options.output)

    try:
        macros = load_macros(
            definition_paths,
            reject_duplicates=options.reject_duplicates,
        )
        expander = MacroExpander(macros, max_depth=options.max_depth)
        expand_file(input_path, output_path, expander)
    except MacroError as exc:
        print(f"macro_expand.py: error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
