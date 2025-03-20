MODULE CRTM_VLIDORT_ChannelInfo
  
  IMPLICIT NONE

  TYPE ChannelInfo_Define
    LOGICAL :: is_wavenumber = .true.
    REAL(8) :: start_wv, end_wv
    REAL(8) :: resolution
    REAL(8), allocatable  :: wavenumber(:)
    REAL(8), allocatable  :: wavelength(:)
    INTEGER :: N_Channels = -1
 
  END TYPE ChannelInfo_Define

CONTAINS

   SUBROUTINE set_usr_wv(ChannelInfo)
      TYPE(ChannelInfo_Define), INTENT(IN OUT) :: ChannelInfo
      REAL(8) :: wl_usr
      INTEGER :: i
   
      ChannelInfo%N_Channels = (ChannelInfo%end_wv-ChannelInfo%start_wv) / ChannelInfo%resolution + 1

      IF (allocated(ChannelInfo%wavenumber) ) deallocate(ChannelInfo%wavenumber)
      IF (allocated(ChannelInfo%wavelength) ) deallocate(ChannelInfo%wavelength)
      allocate(ChannelInfo%wavenumber(1:ChannelInfo%N_Channels) )
      allocate(ChannelInfo%wavelength(1:ChannelInfo%N_Channels) )
      DO i = 1, ChannelInfo%N_Channels
         wl_usr  = ChannelInfo%start_wv + (i-1)*ChannelInfo%resolution
         IF(ChannelInfo%is_wavenumber) THEN
            ChannelInfo%wavenumber(i) = wl_usr
            ChannelInfo%wavelength(ChannelInfo%N_Channels-i+1) =  1.D7/wl_usr
 
         ELSE
           ChannelInfo%wavelength(i) = wl_usr
           ChannelInfo%wavenumber(ChannelInfo%N_Channels-i+1) =  1.D7/wl_usr
         ENDIF
      ENDDO

  End SUBROUTINE set_usr_wv



END MODULE CRTM_VLIDORT_ChannelInfo 
