         !$acc parallel present(mesh) vector_length(128) num_gangs(9750) async(1)
         !$acc loop gang
         do eID = 1 , size(mesh % elements)
            Nx = mesh % elementNxyzFlat(1,eID)
            Ny = mesh % elementNxyzFlat(2,eID)
            Nz = mesh % elementNxyzFlat(3,eID)
            cfStride = NCONS * (Nx + 1) * (Ny + 1) * (Nz + 1)

            !$acc loop vector collapse(3) private(Flux, node, qBase, gradBase, muBase, cfBaseX, cfBaseY, cfBaseZ)
            do k = 0, Nz
               do j = 0, Ny
                  do i = 0, Nx
                     node = 1 + i + (Nx + 1) * (j + (Ny + 1) * k)
                     qBase = mesh % storage % elementNodeOffset(eID) + (node - 1) * NCONS
                     gradBase = mesh % storage % elementGradOffset(eID) + (node - 1) * NGRAD
                     muBase = mesh % storage % elementMuOffset(eID) + (node - 1) * 3

                     call ViscousFlux_STATE(NCONS, NGRAD,                                                &
                                            mesh % storage % QNS(qBase:qBase + NCONS - 1),             &
                                            mesh % storage % U_xNS_pool(gradBase:gradBase + NGRAD - 1), &
                                            mesh % storage % U_yNS_pool(gradBase:gradBase + NGRAD - 1), &
                                            mesh % storage % U_zNS_pool(gradBase:gradBase + NGRAD - 1), &
                                            mesh % storage % mu_NS_pool(muBase), 0.0_RP,               &
                                            mesh % storage % mu_NS_pool(muBase + 1), Flux)

                     cfBaseX = mesh % storage % elementContravariantOffset(eID) + (node - 1) * NCONS
                     cfBaseY = mesh % storage % elementContravariantOffset(eID) + cfStride + (node - 1) * NCONS
                     cfBaseZ = mesh % storage % elementContravariantOffset(eID) + 2 * cfStride + (node - 1) * NCONS

                     !$acc loop seq
                     do eq = 1, NCONS
                        mesh % storage % contravariantFlux_pool(cfBaseX + eq - 1) =                       &
                           - Flux(eq,IX) * mesh % elementJGradXiPool(IX,node,eID)                         &
                           - Flux(eq,IY) * mesh % elementJGradXiPool(IY,node,eID)                         &
                           - Flux(eq,IZ) * mesh % elementJGradXiPool(IZ,node,eID)

                        mesh % storage % contravariantFlux_pool(cfBaseY + eq - 1) =                       &
                           - Flux(eq,IX) * mesh % elementJGradEtaPool(IX,node,eID)                        &
                           - Flux(eq,IY) * mesh % elementJGradEtaPool(IY,node,eID)                        &
                           - Flux(eq,IZ) * mesh % elementJGradEtaPool(IZ,node,eID)

                        mesh % storage % contravariantFlux_pool(cfBaseZ + eq - 1) =                       &
                           - Flux(eq,IX) * mesh % elementJGradZetaPool(IX,node,eID)                       &
                           - Flux(eq,IY) * mesh % elementJGradZetaPool(IY,node,eID)                       &
                           - Flux(eq,IZ) * mesh % elementJGradZetaPool(IZ,node,eID)

                        mesh % storage % QdotNS(qBase + eq - 1) = 0.0_RP
                     end do
                  end do
               end do
            end do

            !$acc loop vector collapse(3) private(node, qdotBase, nodeX, nodeY, nodeZ, cfBaseX, cfBaseY, cfBaseZ)
            do k = 0, Nz
               do j = 0, Ny
                  do i = 0, Nx
                     node = 1 + i + (Nx + 1) * (j + (Ny + 1) * k)
                     qdotBase = mesh % storage % elementNodeOffset(eID) + (node - 1) * NCONS
                     !$acc loop seq
                     do eq = 1, NCONS
                        !$acc loop seq
                        do l = 0, Nx
                           nodeX = 1 + l + (Nx + 1) * (j + (Ny + 1) * k)
                           nodeY = 1 + i + (Nx + 1) * (l + (Ny + 1) * k)
                           nodeZ = 1 + i + (Nx + 1) * (j + (Ny + 1) * l)
                           cfBaseX = mesh % storage % elementContravariantOffset(eID) + (nodeX - 1) * NCONS
                           cfBaseY = mesh % storage % elementContravariantOffset(eID) + cfStride + (nodeY - 1) * NCONS
                           cfBaseZ = mesh % storage % elementContravariantOffset(eID) + 2 * cfStride + (nodeZ - 1) * NCONS

                           mesh % storage % QdotNS(qdotBase + eq - 1) = mesh % storage % QdotNS(qdotBase + eq - 1) + &
                              mesh % nodalHatDPool(Nx,i,l) * mesh % storage % contravariantFlux_pool(cfBaseX + eq - 1) + &
                              mesh % nodalHatDPool(Ny,j,l) * mesh % storage % contravariantFlux_pool(cfBaseY + eq - 1) + &
                              mesh % nodalHatDPool(Nz,k,l) * mesh % storage % contravariantFlux_pool(cfBaseZ + eq - 1)
                        end do
                     end do
                  end do
               end do
            end do
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
            !$acc loop vector collapse(3) private(Flux, node, qBase, nodeX, nodeY, nodeZ, qBaseXiNode, qBaseEtaNode, qBaseZetaNode)
            do k = 0, Nz
               do j = 0, Ny
                  do i = 0, Nx
                     node = 1 + i + (Nx + 1) * (j + (Ny + 1) * k)
                     qBase = mesh % storage % elementNodeOffset(eID) + (node - 1) * NCONS
                     !$acc loop seq
                     do l = 0, Nx
                        nodeX = 1 + l + (Nx + 1) * (j + (Ny + 1) * k)
                        nodeY = 1 + i + (Nx + 1) * (l + (Ny + 1) * k)
                        nodeZ = 1 + i + (Nx + 1) * (j + (Ny + 1) * l)
                        qBaseXiNode = mesh % storage % elementNodeOffset(eID) + (nodeX - 1) * NCONS
                        qBaseEtaNode = mesh % storage % elementNodeOffset(eID) + (nodeY - 1) * NCONS
                        qBaseZetaNode = mesh % storage % elementNodeOffset(eID) + (nodeZ - 1) * NCONS

                        call TwoPointFlux_Selector(mesh % storage % QNS(qBase:qBase + NCONS - 1),              &
                                                   mesh % storage % QNS(qBaseXiNode:qBaseXiNode + NCONS - 1), &
                                                   mesh % elementJGradXiPool(:,node,eID),                     &
                                                   mesh % elementJGradXiPool(:,nodeX,eID), Flux(:,IX))
                        call TwoPointFlux_Selector(mesh % storage % QNS(qBase:qBase + NCONS - 1),                &
                                                   mesh % storage % QNS(qBaseEtaNode:qBaseEtaNode + NCONS - 1), &
                                                   mesh % elementJGradEtaPool(:,node,eID),                      &
                                                   mesh % elementJGradEtaPool(:,nodeY,eID), Flux(:,IY))
                        call TwoPointFlux_Selector(mesh % storage % QNS(qBase:qBase + NCONS - 1),                  &
                                                   mesh % storage % QNS(qBaseZetaNode:qBaseZetaNode + NCONS - 1), &
                                                   mesh % elementJGradZetaPool(:,node,eID),                       &
                                                   mesh % elementJGradZetaPool(:,nodeZ,eID), Flux(:,IZ))

                        !$acc loop seq
                        do eq = 1, NCONS
                           mesh % storage % QdotNS(qBase + eq - 1) = mesh % storage % QdotNS(qBase + eq - 1) &
                              - mesh % nodalSharpDPool(Nx,i,l) * Flux(eq,IX)                                  &
                              - mesh % nodalSharpDPool(Ny,j,l) * Flux(eq,IY)                                  &
                              - mesh % nodalSharpDPool(Nz,k,l) * Flux(eq,IZ)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$acc end parallel loop
