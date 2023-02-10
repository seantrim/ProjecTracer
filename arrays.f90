module arrays
implicit none
real*8, allocatable :: Matrix(:,:),Bvec(:),Row(:),Col(:),MatrixB(:,:),RowB(:),ColB(:)  !!p,q indexing
integer*4, allocatable :: IPIV(:),IPIVB(:),ipvec(:),kpvec(:),parray(:,:),Ttype(:)
integer*4, allocatable :: nempty(:,:),count_type(:,:,:)
real*8, allocatable :: nempty_real(:,:),count_type_real(:,:,:)
real*8, allocatable :: xg(:),zg(:),T(:,:),u(:,:),w(:,:),SF(:,:),vis(:,:) !!i,k indexing
real*8, allocatable :: xr(:),xrr(:),xrrr(:),xrrrr(:),zs(:),zss(:),zsss(:),zssss(:)
real*8, allocatable :: D1(:),D2(:),D3(:),D4(:),Drs(:,:),Drss(:,:),Drrs(:,:),Drrss(:,:),dt_array(:,:),dts_array(:,:)
real*8, allocatable :: factorial(:),D(:,:),dr_power(:),ds_power(:)
real*8, allocatable :: D1_so(:),D2_so(:),Drs_so(:,:)
real*8, allocatable :: rtr(:),str(:),Ttr(:),Ctr(:),DERr_u(:,:,:),DERr_w(:,:,:),DERxr(:,:),DERzs(:,:)
real*8, allocatable :: rtr0(:),rtr1(:),rtr2(:),str0(:),str1(:),str2(:),Cvis(:,:,:),Cbuoy(:,:,:),Tvis(:,:),Tbuoy(:,:)
real*8, allocatable :: Tratio(:,:),Ttr0(:),Ttr1(:),Ttr2(:),DERr_tracer(:,:,:),tracer_space_array(:,:),Textend(:,:)
!real*8, allocatable :: DERr_Textend(:,:,:),Cnew(:,:,:),RaC(:),visC(:),mass(:) !!moved RaC(:) and visC(:) to basics.f90
real*8, allocatable :: DERr_Textend(:,:,:),Cnew(:,:,:),mass(:)
real*8, allocatable :: dt_tr(:),dt_stream(:,:,:),dt_vis(:,:,:),dt_fine(:,:),dt_stream_smoother(:,:,:),dt_stream_diff(:,:,:)
!real*8, allocatable :: conduct(:,:),conduct_r(:,:),conduct_s(:,:),conduct_factor(:),Htype(:) !!conduct_factor and Htype moved to basics.f90
real*8, allocatable :: conduct(:,:),conduct_r(:,:),conduct_s(:,:)
real*8, allocatable :: strain(:,:,:),RHS(:,:,:),residual(:,:,:),error(:,:,:)
real*8, allocatable :: vis_x(:,:,:),vis_z(:,:,:),vis_xx(:,:,:),vis_zz(:,:,:),vis_xz(:,:,:),Matrix_coarse(:,:),B_coarse(:)
real*8, allocatable :: vis_xxpzz(:,:,:),vis_xxmzz(:,:,:),vis_smooth(:,:,:)
real*8, allocatable :: vis_grid(:,:,:),vis_gridf(:,:,:),vis_grid_static(:,:,:),B_coarse_save(:)
real*8, allocatable :: vis_xf(:,:,:),vis_zf(:,:,:),vis_xxf(:,:,:),vis_zzf(:,:,:),vis_xzf(:,:,:),visNN_static(:,:,:)
real*8, allocatable :: vis_xxpzzf(:,:,:),vis_xxmzzf(:,:,:),vis_smoothf(:,:,:)
real*8, allocatable :: ratio(:,:,:),stream_restrict(:,:,:),stream_smooth(:,:,:),yield_strain(:,:,:)
real*8, allocatable :: u_grid(:,:,:),w_grid(:,:,:)
real*8, allocatable :: SFxx(:,:,:),SFzz(:,:,:),SFxxx(:,:,:),SFzzz(:,:,:),SFxxxx(:,:,:),SFzzzz(:,:,:),SFxz(:,:,:),SFxxzz(:,:,:)&
&                      ,SFxxz(:,:,:),SFxzz(:,:,:)
real*8, allocatable :: SI1(:,:,:),SI2(:,:,:),SI4(:,:,:),SI5(:,:,:),SI6(:,:,:)
real*8, allocatable :: SR1(:,:,:),SR2(:,:,:),SR4(:,:,:),SR5(:,:,:),SR6(:,:,:)
real*8, allocatable :: fake(:,:)
real*8, allocatable :: SM1(:,:),SM2(:,:),SM4(:,:),SM5(:,:),SM6(:,:)
real*8, allocatable :: TCheb(:),TCheb_x(:),TCheb_xx(:),bCheb(:),aCheb(:)
real*8, allocatable :: mewCheb(:),mewtCheb(:),vCheb(:),gtCheb(:),cCheb(:)
real*8, allocatable :: TChebS(:),TCheb_xS(:),bChebS(:)
real*8, allocatable :: mewChebS(:),mewtChebS(:),vChebS(:),cChebS(:)
end module arrays
