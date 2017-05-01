function [ ldsDATA ] = setLdsData( Nx, Ny, Nxy, ptLds, ptLoc, p0 )
%setLdsData Sets the external panel loads data

%  ldsDATA = returned data structure

ldsDATA.Nx=Nx;
ldsDATA.Ny=Ny;
ldsDATA.Nxy=Nxy;

ldsDATA.ptLds=ptLds;
ldsDATA.ptLoc=ptLoc;

ldsDATA.p0=p0;