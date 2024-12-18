README for NOAA team

1. The directories (vlidort_NOAA_V1) contain:
   1) PCRTM_VLIDORT_Main
      This is the main folder for driver: pcrtm_omps_simulator.f90
      and modules: GAS_OPT.f90, PCRTM_SW_Aux.f90
      
      
   2) VLIDORT source codes:
      /util 
      /vlidort_def 
      /vlidort_focode 
      /vlidort_main
      
   3) vliodrt testing folders:
      /vlidort_s_test
      /vlidort_sc_test
      /vlidort_v_test
   
    4)   /PCRTM_VLIDORT_Config
          user_control.cfg, this is main control file
          user_wav_inputs.dat : setting wavelength (given a range or a set of wavelength numbers)
          map_pres_grid.dat, 50 level pressure grid
          std_pres_grid.dat, AFGL standard pressure grid
          user_pres_grid.dat_sample, one sample of user pres. grid
          gasTable, spectral band and absorption gases in each band
          lut_PTH2O_layerAvg.dat, gas LUT's layer average values for pres. temp. and H2O vmr
          LUT_P.dat, gas LUT's pressure grid

     5) Input data_and_cotrol
          
       /GAS_XS, trace gas cross section LUT, covering the spectral range from 23700~43000 cm-1 (232~421nm)
       
       /Profile_afgl, AFGL 1~6 atm. profiles
       
        Inputs_User_profile_Geo.dat_NP_J2_midlat: sample for user provided temperature, water vapor and ozone profiles in 101 levels
        profile_201702_290.500.dat.mid.land:  merra2 mid lat. atm/gas profile, providing the default non-user provided gases profiles
        surface_reflectance_sub.bin, 10 types of surface reflectivity file 

     6) Outputs:  
        There are two options: (a) direct upward radiance (normialized); (b) 3 components, which is a more efficients way to construct 
        omps spectrial given different surface reflectnace 
             
2. Compile and running
   1) make clean
   2) ./pcrtm_vlidort_compile.bash v ifort
       change the compiler accordingly, check the makefile in PCRTM_VLIDORT_Main
   3) ./pcrtm_omps_simulator.exe

3. Compile vlidort library 
   1) make -f makefile_vlidort clean
   2) make -f makefile_vlidort
   3) make -f makefile_vlidort install
   4) make -f makefile_vlidort test

4. Compile and run vlidort-crtm  
   1) cd  vlidort-crtm
   2) modify the makefile with the proper CRTM library and the vlidort library pathes
   3) make 
   
