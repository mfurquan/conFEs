#!/bin/bash
gfortran -o test_confess control.f90 GausLeg_table.f90 LinAlg.f90 InForm.f90 test_confess.f90
