<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# GitHub CI/CD scripts for StormRuler🦜
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CI: `ci-ubuntu.yml`
<!----------------------------------------------------------------------------->
Build and run tests on Ubuntu.

[![Ubuntu](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-ubuntu.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-ubuntu.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CI: `ci-macos.yml`
<!----------------------------------------------------------------------------->
Build and run tests on macOS.

[![macOS](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-macos.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-macos.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CI: `ci-windows.yml`
<!----------------------------------------------------------------------------->
Build and run tests on Windows.

[![Windows](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-windows.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-windows.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CI: `ci-pages.yml`
<!----------------------------------------------------------------------------->
Build publish documentation to GitHub pages.

[![Pages](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-pages.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/ci-pages.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CD: `analysis-codeql.yml`
<!----------------------------------------------------------------------------->
Build and analyze with **CodeQL**. *This pipeline does not work as expected due
to the fact that CodeQL uses `clang++` internally to analyze the code,
and Storm fails to compile/be parsed under `clang++`.*

[![CodeQL](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-codeql.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-codeql.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CD: `analysis-coverity.yml`
<!----------------------------------------------------------------------------->
Build and analyze with **Coverity**. *Coverity Build Tool **2022.6.0** is the
current avaliable version, it does not work with Storm properly. Although,
Windows version works fine, so we are using it. At least for now.*

[![Coverity](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-coverity.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-coverity.yml)

--------------------------------------------------------------------------------

<!----------------------------------------------------------------------------->
## CD: `analysis-sonar.yml`
<!----------------------------------------------------------------------------->

Build, run tests, collect coverage, analyze with **SonarCloud** and upload the
coverage data to **Codecov**.

[![SonarCloud](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-sonar.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-sonar.yml)

--------------------------------------------------------------------------------
