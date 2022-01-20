#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using Base.Filesystem
using ArgParse

# copies a GCM environment file to a new directory,
# but 
