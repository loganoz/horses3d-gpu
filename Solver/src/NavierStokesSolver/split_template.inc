         !$acc parallel present(mesh) vector_length(128) num_gangs(9750) async(1)
         !$acc loop gang
         do eID = 1 , size(mesh % elements)
         
            !$acc loop vector collapse(3) private(Flux)
            do k = 0, mesh % elements(eID) % Nxyz(3)  
               do j = 0, mesh % elements(eID) % Nxyz(2)  
                  do i = 0, mesh % elements(eID) % Nxyz(1)

                  call ViscousFlux_STATE( NCONS, NGRAD, mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % U_x(:,i,j,k), & 
                                          mesh % elements(eID) % storage % U_y(:,i,j,k) , mesh % elements(eID) % storage % U_z(:,i,j,k), &
                                          mesh % elements(eID) % storage % mu_ns(1,i,j,k), 0.0_RP, &
                                          mesh % elements(eID) % storage % mu_ns(2,i,j,k), Flux)
                  !$acc loop seq
                  do eq = 1, NCONS
              
                     mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IX)  = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradXi(IX,i,j,k)  &
                                                                                        - Flux(eq,IY) * mesh % elements(eID) % geom % jGradXi(IY,i,j,k)  &
                                                                                        - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)

                     mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IY) = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  &
                                                                                       - Flux(eq,IY) * mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  &
                                                                                       - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)
                  
                     mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IZ) = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradZeta(IX,i,j,k)  &
                                                                                       - Flux(eq,IY) * mesh % elements(eID) % geom % jGradZeta(IY,i,j,k)  &
                                                                                       - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k)
                     ! initialize to 0 to accumulate
                     mesh % elements(eID) % storage % Qdot(eq,i,j,k)  = 0.0_RP
                  end do
            end do               ; end do                ; end do

         call ScalarWeakIntegrals_StdVolumeGreen( mesh % elements(eID) % Nxyz, NCONS, mesh % elements(eID) % storage % contravariantFlux, &
                                                  mesh % elements(eID) % storage % QDot)
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
            !$acc loop vector collapse(3) private(Flux)
            do k = 0, mesh % elements(eID) % Nxyz(3)  
               do j = 0, mesh % elements(eID) % Nxyz(2)  
                  do i = 0, mesh % elements(eID) % Nxyz(1)
                     if ( mesh % elements(eID) % Nxyz(1) > 0 ) then
                        !$acc loop seq
                        do l = 0, mesh % elements(eID) % Nxyz(1)
                           call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,l,j,k), mesh % elements(eID) % geom % jGradXi(:,i,j,k), mesh % elements(eID) % geom % jGradXi(:,l,j,k), Flux(:,IX))

                           !$acc loop seq
                           do eq = 1, NCONS
                              mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) &
                                                                              - NodalStorage(mesh % elements(eID) % Nxyz(1)) % sharpD(i,l) * Flux(eq,IX)
                           end do
                        end do
                     end if

                     if ( mesh % elements(eID) % Nxyz(2) > 0 ) then
                        !$acc loop seq
                        do l = 0, mesh % elements(eID) % Nxyz(2)
                           call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,i,l,k), mesh % elements(eID) % geom % jGradEta(:,i,j,k), mesh % elements(eID) % geom % jGradEta(:,i,l,k), Flux(:,IY))

                           !$acc loop seq
                           do eq = 1, NCONS
                              mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) &
                                                                              - NodalStorage(mesh % elements(eID) % Nxyz(2)) % sharpD(j,l) * Flux(eq,IY)
                           end do
                        end do
                     end if

                     if ( mesh % elements(eID) % Nxyz(3) > 0 ) then
                        !$acc loop seq
                        do l = 0, mesh % elements(eID) % Nxyz(3)
                           call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,i,j,l), mesh % elements(eID) % geom % jGradZeta(:,i,j,k), mesh % elements(eID) % geom % jGradZeta(:,i,j,l), Flux(:,IZ))

                           !$acc loop seq
                           do eq = 1, NCONS
                              mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) &
                                                                              - NodalStorage(mesh % elements(eID) % Nxyz(3)) % sharpD(k,l) * Flux(eq,IZ)
                           end do
                        end do
                     end if

            end do               ; end do                ; end do

         enddo
         !$acc end parallel loop
