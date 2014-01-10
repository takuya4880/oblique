module defstruct
    implicit none
    integer,parameter :: nx=50
    integer,parameter :: nz=100
    integer,parameter :: cox=8
    integer,parameter :: coz=8
    integer,parameter :: nnx=nx*cox
    integer,parameter :: nnz=nz*coz
    integer,parameter :: marg=5
    integer,parameter :: ix=nx+2*marg
    integer,parameter :: iz=nz+2*marg
    integer,parameter :: iix=nnx+2*marg
    integer,parameter :: iiz=nnz+2*marg
    

    type constants 
        integer nx, nz, ix, iz, marg
        integer imx,imz
        double precision dx, dz, dt, wid, hig
        double precision gam, q, a
        double precision gx, gy, gz 
    end type

    type output
        integer mf_params, mf_t, mf_ro, mf_pr, mf_vx
        integer mf_vy, mf_bx, mf_by, mf_az, mf_x, mf_y
        integer mfi_t, mfi_ro, mfi_pr, mfi_vx
        integer mfi_vy, mfi_bx, mfi_by, mfi_az
    end type
    
    type cell
        type(constants) con
        type(output) op
        double precision x(ix), z(iz)
        double precision ro(ix,iz), rovx(ix,iz), rovy(ix,iz), rovz(ix,iz)
        double precision bx(ix,iz), by(ix,iz), bz(ix,iz), e(ix,iz) 
        double precision pr(ix,iz)
        double precision bpot(ix,iz)       
    end type    


    
end module

