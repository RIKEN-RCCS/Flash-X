# Flash-X, a Multiphysics Scientific Software System

Flash-X is a highly composable multiphysics software system that can
be used to simulate physical phenomena in several scientific
domains. It is derived from FLASH, which has a history of being a
community code for several communities. The Flash-X architecture has
been redesigned to be compatible with increasingly heterogeneous
hardware platforms. Part of the redesign is a newly designed
performance portability layer that is language agnostic.

## Libraries

The code of the PARAMESH library for Flash-X is provided as a git submodule by a separate code
repository named PARAMESH-for-Flash-X.  After successful initialization of the submodule
--- see further below ---
the code will appear as a subtree under `Grid/GridMain/AMR` in the `PM4_package` directory.
The Flash-X `setup` and build mechanisms will take care of compiling this code (if needed
for a simulation configuration); PARAMESH code is not organized or built as a separate library.

Some applications and tests use external libraries that are expected to be already installed on the
system where Flash-X is being built and run. The directory locations of such library installations
should be made known to the Flash-X build system by a site-specific (or, as a fallback, OS-specific)
Makefile.h file. See the subdirectories under sites/ .

This applies in particular to the AMReX library. Separate library instances for 1D, 2D, and 3D
should be installed, and the appropriate locations mentioned in Makefile.h . The Flash-X code in
this branch is compatible with AMReX version 23.11 (and with earlier unreleased versions of AMReX that
have f817d776f544410c315cf6007564f88f6b16fc2b and d77d93b060ae4bc7e1e65ac1465ad7fc06a5d4dc in their commit history).
Simulations that use AMReX may not be compatible with versions of AMReX earlier than that.

On the other hand, some applications and tests use INTERNAL libraries. Such libraries are built, as part of the Flash-X setup process, from source code that is located in subdirectories under lib/ . There are two cases for such libraries:

1. **Library source code is included as part of the Flash-X git respository.**

   An example is the sqrt3 library, whose source code is included in lib/sqrt/ .

2. **Library source code must be retrieved from a separate repository.**

   Examples are the THORNADO and WEAKLIB libraries.
   Follow the instructions on submodules to automatically put the source code for these two
   in the right places in subdirectories under lib/.

## Git with Submodules

To prepare for building simulations that use libraries whose code must be retrieved
from a separate git repository, the following modified `git` commands can be used:

- `git pull --recurse-submodules=yes` (in place of the usual `git pull`)
- `git submodule update --init` (additionally, after `git pull`)
- `git submodule update --init source/Grid/GridMain/AMR/Paramesh4/PM4_package`
  (variant of the previous item, in case you want to only get the Paramesh package submodule)

## Git/Testing Workflow

The current rules for collaborating via GitHub are as follows:

Contributors with
read only permission to the Flash-X code repository should use the following
guidelines to create a pull request:

1. Create a fork.
2. Make your changes.
3. Create a PR to the **staged** branch whenever you wish.
   Give your PR a title that begins with the word "DRAFT".
   This will allow any discussion about the pull
   request to be conducted on github.
4. When you are ready for the pull request to be accepted, merge from **main**
   into your forked code, to ensure that your fork is not out of sync.
4. If a merge conflict occurs when merging **main** into the feature branch,
   _do not_ attempt to resolve conflicts using the  GitHub web interface - such an attempt can results in an unintended merge to **main**.
5. Run a local version of your test suite and make sure everything
   passes.
6. Make sure your latest commit has been pushed.
7. Remove "DRAFT" from your pull request name. If no further problems
   are found, this will cause the PR
   to be merged. The test suite is run at night if one of more
   PRs have been merged into the **staged** branch. PRs that come in
   before 6PM CST are more likely to be included in that night's test.
   Monitor the repo to see whether your PR was merged and the test suite passed.
   A comment will be added to your PR if the test suite failed.
8. If the test suite passes, a composite PR will be created from
   **staged** into **main**, and you won't have to do anything more. This will
   likely happen the day the test suite passed.
9. If the test suite fails, it is expected that you will prioritize resolving the
   failure. Note that the merged and colliding code will be available in the staged branch.
   You can copy that code into a local working copy to resolve the issue. **Please note that you
   should never make any commit into the staged branch.**
   If the test suite passes, you can reissue a PR and ask for a test suite run by leaving a comment in the PR.
   If the failures continue, we abandon the stage branch at the end of the day, and PRs have to be created again.
   If we determine that the interoperability is compromised, someone from the core team might have
   to get involved to help resolve.

Contributors with write permission should create a feature branch from the main branch
instead of a fork. The remainder of the workflow remains the same.

## Code Formatting

Currently Flash-X source code uses different formatting styles for different code units:

1. **fprettify** (https://github.com/pseewald/fprettify): With 3 white space indentation style.
2. **fxprettify** (tools/fileTools/fxprettify.py): A wrapper around fprettify that implements emacs style formatting. 

Until strict formatting rules are established, we recommend using a style that produces minimum whitespace 
differences when editing/adding code to an existing unit. For new units we recommend using **fprettify** with
the following setting,

```
fprettify -i3 -w1 *.F90
```

## Special Syntax Highlighting

To enable syntax for **Flash-X** specific keywords implement following settings based on the editor:

### **VIM**

- Create `$HOME/.vim/after/syntax/fortran.vim` and add:

```
" Custom FORTRAN keywords
syn match customKeyword '\<*\zs@MACRO\>'
syn match customKeyword '\<*\zs@M\>'
syn match customKeyword '\<*\zs!!REORDER\>'
syn match customKeyword '\<*\zs!!NOVARIANTS\>'

" Definitions
hi def link customKeyword Keyword
```

- Create `$HOME/.vimrc` and add:
```
" Turn on syntax
syntax on

" Set file type
autocmd BufNewFile,BufRead *.F90-mc set filetype=fortran
```

## Tests
Test specifications for individual simulations are included under ``*/tests/tests.yaml`` files in each simulation directory under ``source/Simulation/SimulationMain``. New tests should be added as enteries in the prescribed YAML format before including it as a part of suites on different platforms.

Testing and maintainence of the code is implemented using command line tools available in Flash-X-Test repository: https://github.com/Flash-X/Flash-X-Test

Please refer to the instructions there to setup your own testing infrastructure. Also take a look at ``sites/ganon_jenkins/UnitTests.suite`` for an example of publicly available test suite which can be edited to enable code coverage for new modules.

Testing servers:

- Argonne, GCE:

  FlashTest server for running tests on `staged` branch

  - GCC   - https://jenkins.cels.anl.gov/job/Flash-X-staged_GCC
  - Intel - https://jenkins.cels.anl.gov/job/Flash-X-staged_Intel

  FlashTestView

  - GCC - https://web.cels.anl.gov/projects/FLASH5/testsuite/home.py?target_dir=/nfs/pub_html/gce/projects/FLASH5/output/staged_gcc
  - Intel - https://web.cels.anl.gov/projects/FLASH5/testsuite/home.py?target_dir=/nfs/pub_html/gce/projects/FLASH5/output/staged_intel

- Ganon:

  FlashTest - http://ganon.device.utk.edu:8080 or http://ganon2.device.utk.edu:8080

- Compatibility with AMReX

  Testing with latest AMReX developments happens here - https://jenkins-gce.cels.anl.gov/job/Flash-X-amrex_GCC

## Citation

Please use the following for citing Flash-X

```
@article{Flash-X-SoftwareX,
title = {Flash-X: A multiphysics simulation software instrument},
journal = {SoftwareX},
volume = {19},
pages = {101168},
year = {2022},
issn = {2352-7110},
doi = {https://doi.org/10.1016/j.softx.2022.101168},
url = {https://www.sciencedirect.com/science/article/pii/S2352711022001030},
author = {Anshu Dubey and Klaus Weide and Jared O’Neal and Akash Dhruv and Sean Couch and J. Austin Harris and Tom Klosterman and Rajeev Jain and Johann Rudi and Bronson Messer and Michael Pajkos and Jared Carlson and Ran Chu and Mohamed Wahib and Saurabh Chawdhary and Paul M. Ricker and Dongwook Lee and Katie Antypas and Katherine M. Riley and Christopher Daley and Murali Ganapathy and Francis X. Timmes and Dean M. Townsley and Marcos Vanella and John Bachan and Paul M. Rich and Shravan Kumar and Eirik Endeve and W. Raphael Hix and Anthony Mezzacappa and Thomas Papatheodore},
}
```
