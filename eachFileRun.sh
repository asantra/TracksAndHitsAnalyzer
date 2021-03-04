#! /bin/bash

root -l -b << EOF
  .L process_lxtrees_v3.C++
  process_lxtrees_v3("$1", "$2")
EOF
