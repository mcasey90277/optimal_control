function [q, qd, L, A] = sleigh_lagrangian(p)
%% Purpose:
%
%   The sleigh as mech.lagrange_dalembert wants it: generalized
%   coordinates q = (x_c, y_c, th) of the CENTER OF MASS, the free-body
%   Lagrangian L = m/2 |rd_c|^2 + J/2 thdot^2 (no potential), and the
%   knife-edge constraint at the blade point r_b = r_c - a e1: zero
%   sideways velocity, rd_b . e2 = -xd_c sin th + yd_c cos th - a thdot = 0,
%   so A = [-sin th, cos th, -a]. Nothing is solved here.
%
%% Inputs:
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  q                        [3 x 1] sym             [xc; yc; th]
%
%  qd                       [3 x 1] sym             [xcd; ycd; thdot]
%
%  L                        sym                     Lagrangian
%
%  A                        [1 x 3] sym             Constraint row
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: print the objects for the defaults
     [q, qd, L, A] = sleigh_lagrangian(sleigh_params(1, 0.1, 0.2));
     disp(q.');  disp(qd.');  disp(L);  disp(A);
     return
end
syms xc yc th xcd ycd thdot real            % not `thd`: a Signal Processing Toolbox function shadows that name
q  = [xc; yc; th];
qd = [xcd; ycd; thdot];
L  = p.m/2*(xcd^2 + ycd^2) + p.J/2*thdot^2;
A  = [-sin(th), cos(th), -p.a];
end
