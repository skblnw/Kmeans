program Kmean

use omp_lib

integer,parameter :: K=8, N=230, M=1
double precision :: lamb = 0 !0
double precision,parameter :: a = 0.0001 !0.000025
!double precision :: lamb
double precision :: C(N,N), D(N,N), WCSS2, WCSS3, WCSS, FWCSS(N), AWCSS(N), BWCSS(N), Min_WCSS, Min_WCSS2, Min_WCSS3
character(len=40) :: fname, fnameD

! S(:,k) : indices of atoms in cluster k until value = 0 
! S(:,K) = (/2 5 1  0 0 0 0 0...0 /) means cluster k contains atoms with indices 2 , 5 1,             
!       and    S_size(k) = 3
! S_size(k) = # of atoms in cluster k = # of non-zeros in S(:,k)
! S_label(i) = index of cluster where atom i is .

! C : essential r-covariance matrix 
! D : r-covariance matrix
! WCSS:  within-cluster sum of squares


integer :: S(N,K), S_size(K), S_label(N),T(K)
integer :: S_new(N,K), S_size_new(K), S_label_new(N)
logical :: FINISHED,FINIED


integer ::  atom_i, k_iter,tt, iter, Max_Iter = 100, min_indx(1),indx;     

! coved.dat is the essensial subspace covariance matrix (natm*natm). Calculated by adding up essensial eigenvector*eigenvalue*eigenvector and x+y+z
fname = "coved.dat";
! covd.dat is the (nearly) full covariance matrix (natm*natm). Elements are C(i,j) = x(i)*x(j) + y(i)*y(j) + z(i)*z(j)
fnameD ="covd.dat";
call Read_Data(N,C,fname)
call Read_DataD(N,D,fnameD)

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'KMEANS'
write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  Do PCA-EDCG with an additional geometrical factor.'

write ( *, '(a)' ) ' '
write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

write(*,"(A,I0,A,I0)") "number of atom (N): ", N,  " | number of cluster (K): ", K

FINIED = .FALSE. ; tt=0
DO WHILE ( .NOT. FINIED )
    write(*,"(A, F10.6)"), "lamb=" ,lamb
    tt = tt+1
    print *, tt

    do atom_i=1,N
        !print *, "i=" ,atom_i

        call Cal_initial(atom_i,N,K,C,D,T,S,S_size,S_label)

        WCSS2 = Cal_WCSSed( K,N,C,S,S_size,S_label )
        WCSS3 = Cal_WCSSsp( K,N,D,S,S_size,S_label )
        WCSS = WCSS2 + lamb * wcss3
        !write(*,*) "Init is setup well and saved!" 
        !write(*,*) "  WCSS0 = ", WCSS

        FINISHED = .FALSE. ;  iter = 0;

        DO WHILE ( .NOT. FINISHED )
            iter = iter +1
            call   Update_label(lamb, K,N,C,D,S, S_size, S_label_new)
            call   Update_S(K,N, S_new, S_size_new, S_label_new)
            ! convergence iff  S_label = A_label_new
            if (  all (S_label_new .EQ. S_label ) ) then
                !print *, " Converged";
                FINISHED =.TRUE. 
            endif 

            S_label = S_label_new;
            S = S_new ;
            S_size = S_size_new ;

            WCSS2 = Cal_WCSSed( K,N,C,S,S_size,S_label ) 
            WCSS3 = Cal_WCSSsp( K,N,D,S,S_size,S_label )
            WCSS = WCSS2 + lamb * wcss3
            !print *, "iter=", iter, " WCSS=", WCSS

            if ( iter > Max_Iter  )   exit 
        END DO

        !print *, " final WCSS=" , WCSS

        FWCSS(atom_i)=WCSS !given i and lamb
        AWCSS(atom_i)=WCSS2
        BWCSS(atom_i)=WCSS3
    enddo

    Min_WCSS = MINVAL(FWCSS)
    min_indx = minloc(FWCSS(1:N))
    indx=min_indx(1)
    Min_WCSS2=AWCSS(indx)
    Min_wcss3=Bwcss(indx)
    !print *, " minWCSS=" , Min_WCSS
    call Output_FWCSS(lamb,K,N,FWCSS,Min_WCSS,S, S_size, 'FWCSS.dat');
    call Output_AWCSS(lamb,K,N,FWCSS,Min_WCSS2,AWCSS,S, S_size, 'AWCSS.dat');
    call Output_BWCSS(lamb,K,N,FWCSS,Min_WCSS3,BWCSS,S, S_size, 'BWCSS.dat');
    call Output_F(lamb,K,N,FWCSS,Min_WCSS,S, S_size, 'F.dat');
    call Output_A(lamb,K,N,FWCSS,Min_WCSS2,AWCSS,S, S_size, 'A.dat');
    call Output_B(lamb,K,N,FWCSS,Min_WCSS3,BWCSS,S, S_size, 'B.dat');

    atom_i=indx
    call Cal_initial(atom_i,N,K,C,D,T,S,S_size,S_label)
    WCSS2 = Cal_WCSSed( K,N,C,S,S_size,S_label ) 
    WCSS3 = Cal_WCSSsp( K,N,D,S,S_size,S_label )
    WCSS = WCSS2 + lamb * wcss3
    call Output_Data( lamb,atom_i,K,N,S, S_size,  'init.dat' );

    FINISHED = .FALSE. ; iter = 0; 

    DO WHILE ( .NOT. FINISHED ) 
        iter = iter +1 
        call   Update_label(lamb, K,N,C,D,S, S_size, S_label_new)
        call   Update_S(K,N, S_new, S_size_new, S_label_new)
        ! convergence iff  S_label = A_label_new
        if (  all (S_label_new .EQ. S_label ) ) then 
            !print *, " Converged";
            FINISHED =.TRUE. 
        endif 

        S_label = S_label_new; 
        S = S_new ;
        S_size = S_size_new ; 

        WCSS2 = Cal_WCSSed( K,N,C,S,S_size,S_label ) 
        WCSS3 = Cal_WCSSsp( K,N,D,S,S_size,S_label )
        WCSS = WCSS2 + lamb * wcss3
        !print *, "iter=", iter, " WCSS=", WCSS

        if ( iter > Max_Iter  )   exit 
    END DO

    call Output_Data( lamb,atom_i, K,N,S, S_size,  'site.dat' );
    lamb = lamb + a 

    if ( tt > M )   then
        print *, " DONE";
        FINIED =.TRUE. 
    endif 

END DO


contains


subroutine Read_Data(N,C,fname)
    integer , intent(in) :: N            
    character(len=*) :: fname
    double precision :: C(N,N)
    ! FANJUN .........
    double precision :: coved(N,N)
    integer :: ii, jj

    open(11,file=fname)
    do jj=1,N
      do ii=jj,N 
        read(11,*) coved(ii,jj)
        C(jj,ii) = coved(ii,jj)
        C(ii,jj) = C(jj,ii) 
      end do
    end do 
    close(11)
end subroutine  Read_Data

subroutine Read_DataD(N,D,fnameD)
    integer , intent(in) :: N            
    character(len=*) :: fnameD
    double precision :: D(N,N)
    ! FANJUN .........
    double precision :: covedd(N,N)
    integer :: ii, jj

    open(11,file=fnameD)
    do jj=1,N
      do ii=jj,N 
        read(11,*) covedd(ii,jj)
        D(jj,ii) = covedd(ii,jj)
        D(ii,jj) = D(jj,ii) 
      end do
    end do 
    close(11)
end subroutine  Read_DataD


subroutine Output_Data(lamb,atom_i,K,N,S, S_size, fname)
    integer , intent(in) :: N,K,atom_i
    double precision :: lamb 
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K) 

    open(11,file=fname,form="formatted",access="append" )
    write(11,"(A, F10.6)") "lambda=", lamb
    write(11,"(A, I0)") "i=", atom_i
    write(11,"(A, F10.6)") "WCSS=", WCSS
    do k_iter = 1, K
        write(11,"(A, I0, A)",advance="no") "S", k_iter, ": "
        do i_iter = 1, S_size(k_iter)
            write(11,"(I4)",advance="no") S( i_iter, k_iter )
        enddo
        write(11,*) ""
    enddo
    write(11,*) ""
    close(11)
end subroutine  Output_Data

subroutine Output_FWCSS(lamb,K,N,FWCSS,Min_WCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS 
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)
    open(11,file=fname,form="formatted",access="append" )
    write(11,*)  Min_WCSS
    close(11)
end subroutine  Output_FWCSS

subroutine Output_AWCSS(lamb,K,N,FWCSS,Min_WCSS2,AWCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS2, AWCSS(N)
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)

    open(11,file=fname,form="formatted",access="append" )
    write(11,*)  Min_WCSS2
    close(11)
end subroutine  Output_AWCSS
        
subroutine Output_BWCSS(lamb,K,N,FWCSS,Min_WCSS3,BWCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS3, BWCSS(N)
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)

    open(11,file=fname,form="formatted",access="append" )
    write(11,*)  Min_WCSS3
    close(11)
end subroutine  Output_BWCSS

subroutine Output_F(lamb,K,N,FWCSS,Min_WCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS 
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)

    open(11,file=fname,form="formatted",access="append" )
    write(11,*) "lamb=", lamb
    do j=1, N
        write(11,*)  FWCSS(j)
    enddo
    close(11)
end subroutine  Output_F

subroutine Output_A(lamb,K,N,FWCSS,Min_WCSS2,AWCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS2, AWCSS(N)
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)

    open(11,file=fname,form="formatted",access="append" )
    write(11,*) "lamb=", lamb
    do j=1, N
        write(11,*)  AWCSS(j)
    enddo
    close(11)
end subroutine  Output_A
        
subroutine Output_B(lamb,K,N,FWCSS,Min_WCSS3,BWCSS,S, S_size, fname)
    integer , intent(in) :: N,K
    double precision :: lamb,FWCSS(N) ,Min_WCSS3, BWCSS(N)
    integer :: j
    character(len=*) :: fname
    integer, intent(in) :: S(N,K), S_size(K)

    open(11,file=fname,form="formatted",access="append" )
    write(11,*) "lamb=", lamb
    do j=1, N
        write(11,*)  BWCSS(j)
    enddo
    close(11)
end subroutine  Output_B


function Cal_dist(lamb,i,N,K,C,D) result(dist)
    integer :: i,j,N,K,k_iter
    double precision :: C(N,N), D(N,N),lamb,dist(N)
    do j=1,N
        dist(j) = C(i,i) - 2 * C(i,j) + C(j,j) + lamb * D(i,i) - 2 * lamb * D(i,j) + lamb * D(j,j) 
    enddo
end function

subroutine Cal_initial(atom_i,N,K,C,D,T,S,S_size,S_label)
    integer :: S(N,K),S_label(N),S_size(K),T(K)
    double precision :: C(N,N), D(N,N)
    integer :: N,K,atom_i,j,kk,temp
    !integer:: min_indx(N),max_indx(N)
    integer :: max_indx(1),min_indx(1)
    double precision :: dist_min(N), dist_tmp(N), L(K)
    double precision, dimension(K, N) :: dist

    ! do kk=1,K
    !     do j=1,N
    !         dist(kk,j) = 0.0d0
    !     enddo
    ! enddo
    T(1) = atom_i
    do kk=1,K-1
        dist_tmp = Cal_dist(lamb,T(kk),N,K,C,D)
        do j=1,N
            dist(kk,j) = dist_tmp(j)
            dist_min(j) = minval(dist(1:kk,j))
        enddo
        max_indx = maxloc(dist_min(1:N))
        T(kk+1) = max_indx(1)
    enddo

    do j=1,N
        do kk=1,K
            L(kk) = C(T(kk),T(kk)) - 2 * C(T(kk),j) + C(j,j) + lamb * ( D(T(kk),T(kk)) - 2 * D(T(kk),j) + D(j,j) )
        enddo
        min_indx=minloc(L(1:K))
        S_label(j)=min_indx(1)
    enddo

    do K_iter = 1, K
        temp=0
        do j = 1, N
            if (S_label(j) == K_iter) then
                temp = temp + 1;
                S(temp, K_iter) = j;
            endif
        enddo
        S_size(K_iter) = temp;
    enddo
end subroutine

function Cal_distsp( i,k_iter,K,N,D,S,S_Size ) result(dist3)
! calculate SPACIAL distance of atom i to mean of cluster k_iter 
    integer, intent(in) :: i, k_iter, K, N
    double precision, intent(in) :: D(N,N)
    integer, intent(in) :: S(N,K), S_size(K)

    integer :: j, j2
    double precision :: dist3, ss, sss, tmp1, tmp2

    tmp1=0.d0; tmp2=0.d0;
    sss= S_size(k_iter) * S_size(k_iter)
    do j = 1, S_size(k_iter)
        do j2 = 1, S_size(k_iter)
            tmp1= tmp1 + D(S(j, k_iter), S(j2, k_iter) )
        enddo
    enddo
    tmp1 = tmp1 / sss

    ss= S_size(k_iter)
    do j = 1, S_size(k_iter)
        tmp2 = tmp2 -  D(i, S(j,k_iter))
    enddo
    tmp2 = 2 * tmp2 / ss

    dist3 = D(i,i) + tmp2 + tmp1
end function Cal_distsp

function Cal_WCSSsp( K,N,D,S,S_size,S_label ) result (WCSSsp)
    integer, intent(in) :: K, N
    double precision, intent(in) :: D(N,N)
    integer, intent(in) :: S(N,K), S_size(K), S_label(N)

    integer :: i,j,j2, k_iter, min_indx(1), new_site_indx
    double precision ::  ss, WCSSsp, tmp 

    WCSSsp = 0.d0;

    do k_iter = 1 , K  ! visit cluster         
        tmp = 0.d0; 
        do j = 1, S_size(k_iter)  ! visit atom in cluster 
            i = S(j, k_iter) ! atom label 
            tmp = tmp + Cal_distsp(i,k_iter, K, N, D, S, S_Size)           
        enddo 
        WCSSsp = WCSSsp + tmp;
    enddo 
end function Cal_WCSSsp

function Cal_disted( i, k_iter, K, N, C, S, S_Size ) result(dist2)
! calculate distance of atom i to mean of cluster k_iter
    integer, intent(in) :: i, k_iter, K, N
    double precision, intent(in) :: C(N,N)
    integer, intent(in) :: S(N,K), S_size(K)

    integer :: j,j2
    double precision :: dist2,ss,sss,tmp1,tmp2

    tmp1=0.d0; tmp2=0.d0;
    sss= S_size(k_iter) * S_size(k_iter)
    do j = 1, S_size(k_iter)
        do j2 = 1, S_size(k_iter)
            tmp1= tmp1 + C(S(j, k_iter), S(j2, k_iter))
        enddo
    enddo
    tmp1 = tmp1 / sss

    ss = S_size(k_iter)
    do j = 1, S_size(k_iter)
        tmp2 = tmp2 - C(i, S(j,k_iter))
    enddo
    tmp2 = 2 * tmp2 / ss

    dist2 = C(i,i) + tmp2 + tmp1
end function Cal_disted

function Cal_WCSSed( K,N,C,S,S_size,S_label ) result (WCSS2)
    integer, intent(in) :: K, N
    double precision, intent(in) :: C(N,N)
    integer, intent(in) :: S(N,K), S_size(K), S_label(N)

    integer :: i,j,j2, k_iter, min_indx(1), new_site_indx
    double precision :: dist(K), ss, WCSS2, tmp 

    WCSS2 = 0.d0;

    do k_iter = 1 , K  ! visit cluster         
        tmp = 0.d0; 
        do j = 1, S_size(k_iter)  ! visit atom in cluster 
            i = S(j, k_iter) ! atom label                
            tmp = tmp + Cal_disted(i, k_iter, K, N, C, S, S_size )
        enddo
        WCSS2 = WCSS2 + tmp;
    enddo
end function Cal_WCSSed

function Cal_dist2mean(i,k_iter, lamb, K,N,C,D, S, S_Size) result(dist)
! calculate distance of atom i to mean of cluster k_iter 
    integer, intent(in) :: i, k_iter, K, N
    double precision :: lamb 
    double precision, intent(in) :: C(N,N), D(N,N)
    integer, intent(in) :: S(N,K), S_size(K)

    integer :: j,j2
    double precision :: dist, ss,sss,tmp1,tmp2

    tmp1=0.d0; tmp2=0.d0;
    sss= S_size(k_iter) * S_size(k_iter)
    do j = 1, S_size(k_iter) 
        do j2 = 1, S_size(k_iter)
            tmp1= tmp1 +  lamb * D(S(j, k_iter), S(j2, k_iter) )  +  C(S(j, k_iter), S(j2, k_iter) ) 
        enddo  
    enddo 
    tmp1 = tmp1 / sss

    ss= S_size(k_iter) 
    do j = 1, S_size(k_iter) 
        tmp2 = tmp2  -  lamb * D(i, S(j, k_iter) )  -  C(i, S(j,k_iter)) 
    enddo   
    tmp2 = 2 * tmp2 / ss 

    dist = lamb * D(i,i) +  C(i,i) + tmp2 + tmp1 
end function Cal_dist2mean


subroutine  Update_label( lamb, K,N,C,D,S, S_size,S_label_new )
    use omp_lib

    integer, intent(in) :: K, N
    double precision :: lamb 
    double precision, intent(in) :: C(N,N), D(N,N)

    integer, intent(in) :: S(N,K), S_size(K)
    integer :: S_label_new(N)

    integer :: i, k_iter, min_indx(1)
    double precision :: dist(K)

    logical :: debug = .true.
    
    real    ( kind = 8 ) time_begin
    real    ( kind = 8 ) time_elapsed
    real    ( kind = 8 ) time_stop

    WCSS = 0.d0; 

    
    time_begin = omp_get_wtime ( )
    
    !$omp parallel &
    !$omp shared ( N, K, lamb, C, D, S, S_size, min_indx, dist, S_label_new ) &
    !$omp private ( i, k_iter )

    !$omp do
    do i = 1 , N
        ! calculate dist from atom i to all means using old info
        do k_iter = 1, K ! current cluster
            dist(k_iter) = Cal_dist2mean( i,k_iter, lamb, K,N, C,D, S, S_Size)
        enddo
        min_indx = minloc( dist(1:K) ) 
        S_label_new(i) = min_indx(1);
    enddo
    !$omp end do

    !$omp end parallel
    
    time_stop = omp_get_wtime ( )
    
    time_elapsed = time_stop - time_begin
!    write ( *, '(a,g14.6)' ) '  Elapsed time dT = ', time_elapsed
end subroutine   Update_label

subroutine  Update_S(K,N, S_new, S_size_new, S_label_new)
    integer, intent(in) :: K, N
    integer :: S_new(N,K), S_size_new(K), S_label_new(N)
    integer :: i,j, k_iter, temp
    logical :: debug = .true.

    temp=0; j=0;
    
    do K_iter = 1, K
        do i = 1, N
            if (S_label_new(i) == K_iter) then
                temp = temp + 1;
                j = j+1;
                S_new(j, K_iter) = i;
            endif
        enddo
        S_size_new(K_iter) = temp;
        j=0.d0
        temp=0
    enddo
end subroutine   Update_S

end program
