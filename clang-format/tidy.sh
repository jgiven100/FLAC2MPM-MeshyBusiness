#!/bin/bash
find ../ -iname *.h -o -iname *.cpp | xargs clang-format -style='LLVM' -i