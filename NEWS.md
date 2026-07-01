# Changelog

This file documents notable changes in the repository.

## Changes since the merge develop->main on Feb 16, 2026

#### Added
- Added integration with codecov to obtain coverage reports for each PR.
- Added the possibility to write and read partitioning files to avoid the need to 
  partition the mesh at each run.

#### Changed
- The horses3d and horses3d-gpu repositories are now hosted in a GitHub organization
  https://github.com/horses-framework/
- Modified the preprocessing of HOPR HDF5 meshes. The read-in is now MPI-parallel
  when choosing `partitioning = SFC`
