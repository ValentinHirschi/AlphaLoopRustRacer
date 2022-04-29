        subroutine load_inputs(n_inputs, lm, sglm, extm, sgextm, target_res)
            implicit none
            double precision lm(%(n_lm)d,0:3,%(n_inputs)d)
            double precision sglm(%(n_sglm)d,0:3,%(n_inputs)d)
            double precision extm(%(n_extm)d,0:3,%(n_inputs)d)
            double precision sgextm(%(n_sgextm)d,0:3,%(n_inputs)d)
            integer n_lm, n_ext, n_sg_ext, n_inputs
            double precision target_res(%(n_inputs)d)
            integer I,J
            double precision buff(0:(%(n_lm)d+%(n_extm)d+%(n_sgextm)d)*4)

            OPEN(967, FILE='./%(input_file_name)s', ERR=976, STATUS='OLD', ACTION='READ')
                READ(967,*,END=978) n_lm, n_ext, n_sg_ext, n_inputs
                DO I=1,n_inputs
                    READ(967,*,END=978) (buff(J),J=0,(n_lm+n_ext+n_sg_ext)*4)
                    DO J=0,n_lm-1
                        lm(J+1,:,I) = buff(J*4:(J+1)*4-1)
                    ENDDO
                    DO J=0,n_lm+n_ext-1
                        sglm(J+1,:,I) = buff(J*4:(J+1)*4-1)
                    ENDDO
                    DO J=n_lm,n_lm+n_ext-1
                        extm(J-n_lm+1,:,I) = buff(J*4:(J+1)*4-1)
                    ENDDO
                    DO J=n_lm+n_ext,n_lm+n_ext+n_sg_ext-1
                        sgextm(J-n_lm-n_ext+1,:,I) = buff(J*4:(J+1)*4-1)
                    ENDDO
                    target_res(I) = buff((n_lm+n_ext+n_sg_ext)*4)
                ENDDO
                GOTO 978
976         CONTINUE
                STOP 'Could not read the inputs.txt data.'
978         CONTINUE
            CLOSE(967)

        end subroutine load_inputs

        subroutine race(in_lm, in_sglm, in_extm, in_sgextm, out_res)
            implicit none
            double precision in_lm(0:%(n_lm)d-1,0:3)
            double precision in_sglm(0:%(n_sglm)d-1,0:3)
            double precision in_extm(0:%(n_extm)d-1,0:3)
            double precision in_sgextm(0:%(n_sgextm)d-1,0:3)
            double precision out_res

            %(float_type)s lm(0:%(n_lm)d-1,0:3)
            %(float_type)s sglm(0:%(n_sglm)d-1,0:3)
            %(float_type)s extm(0:%(n_extm)d-1,0:3)
            %(float_type)s sgextm(0:%(n_sgextm)d-1,0:3)
            %(float_type)s res

            %(float_type)s em(0:%(n_edges)d-1,0:3)
            %(float_type)s sd(0:%(n_edges)d-1,0:%(n_edges)d-1)
            %(float_type)s em_osE(0:%(n_edges)d-1)
            %(float_type)s num, denom

            double precision pi
            parameter (pi = acos(-1.0d0))
            
            integer I,J
            include 'globals.inc'

            if (DBG) then
                do I=0, %(n_lm)d-1
                    write(*,*) "lm(",I,")=",(in_lm(I,J),J=0,3)
                enddo
                do I=0, %(n_sglm)d-1
                    write(*,*) "sglm(",I,")=",(in_sglm(I,J),J=0,3)
                enddo
                do I=0, %(n_extm)d-1
                    write(*,*) "extm(",I,")=",(in_extm(I,J),J=0,3)
                enddo
                do I=0, %(n_sgextm)d-1
                    write(*,*) "sgextm(",I,")=",(in_sgextm(I,J),J=0,3)
                enddo
            endif

!           Cast all inputs to arbitrary precision
            do J=0,3
                do I=0,%(n_lm)d-1
                    lm(I,J) = REAL(in_lm(I,J), KIND=16)
                enddo
                do I=0,%(n_sglm)d-1
                    sglm(I,J) = REAL(in_sglm(I,J), KIND=16)
                enddo
                do I=0,%(n_extm)d-1
                    extm(I,J) = REAL(in_extm(I,J), KIND=16)
                enddo
                do I=0,%(n_sgextm)d-1
                    sgextm(I,J) = REAL(in_sgextm(I,J), KIND=16)
                enddo
            enddo

            res = 0.0E0_16

!           Evaluae all common quantities
%(warmup_code)s
            
!           Now evaluate all cuts
%(evaluate_ltd_cut)s

!           Finalise computation
            out_res = res
%(wrapup_code_f128)s

        end subroutine race

        program test_cpu_time
            implicit none
            integer i, j
            integer n_inputs
            double precision lm(%(n_lm)d,0:3,%(n_inputs)d)
            double precision sglm(%(n_sglm)d,0:3,%(n_inputs)d)
            double precision extm(%(n_extm)d,0:3,%(n_inputs)d)
            double precision sgextm(%(n_sgextm)d,0:3,%(n_inputs)d)
            double precision res(%(n_inputs)d), target_res(%(n_inputs)d)
            double precision max_diff
            integer max_diff_index
            real :: start, finish
            integer n_runs
            parameter (n_runs =1)

            include 'globals.inc'

            call load_inputs(n_inputs, lm, sglm, extm, sgextm,target_res)

            call cpu_time(start)
            do i=1,n_runs
                do j=1,n_inputs
                    call race(lm(1,0,j),sglm(1,0,j),extm(1,0,j),sgextm(1,0,j),res(j))
                    if (DBG) then
                        write(*,*) "Result :", res(j)
                        write(*,*) "Target :", target_res(j)
                    endif
                enddo
                if (DBG) then
                    stop
                endif
            enddo
            call cpu_time(finish)
            write(*,*) ((finish-start)*1000000.0d0)/(n_runs*n_inputs)
            max_diff = -1.0d0
            do J=1,n_inputs
                if (dabs((res(j)-target_res(j))/target_res(j)).gt.max_diff) then
                    max_diff = dabs((res(j)-target_res(j))/target_res(j))
                    max_diff_index = J
                endif
            enddo
            write(*,*) max_diff*100.0d0
            write(*,*) 'Max diff of ',max_diff*100.0d0,'%% for input #',max_diff_index,' with res=',res(max_diff_index),' and target_res=',target_res(max_diff_index)
        end program test_cpu_time
