function [conn_lines,Connected_Pts,line_num]=Connect_pts(p,q,b,list_ind,conn_lines,Connected_Pts,line_num,Nx,Ny,k,gradDX,gradDY,close_nonadj_pts,count_nonadj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connect_pts function connects a point pair 
% Params:
%   p,q are the point pair that needs to be connected 
%   b is the clumped nuclei boundary point coordinates
%   list_ind is the index number of candidate points
%   conn_lines records the connected lines
%   Connected_Pts records the connected point pairs information
%   line_num is the number of connected point pairs
%   Nx,Ny are the image width, height
%   gradDX,gradDY are the image gradient value along the x-axis and y-axis
%   k is the curvature value of clumped nuclei boundary points
%   curve is the connecting points coordinates 
%   close_nonadj_pts stores those close nonadjacent point pairs
%   count_nonadj is the number of those close nonadjacent point pairs
% Return:
%   conn_lines,Connected_Pts,line_num record the update information
%   =======================================================================================
%   Copyright (C) 2018  Xiaoyuan Guo
%   Email: xiaoyuan.guo@emory.edu
%   =======================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bpt = @(pt) b(list_ind(pt),:);
curve=[];
global bw_poly;

if p==q
    return;
end
[p,q]=Swap_pq(p,q);

if inbackground(bw_poly, bpt(p),bpt(q))
    return;
end

if Intercept_With_Other_Lines(p,q,conn_lines,line_num,b,list_ind)
     return
end

x1 = bpt(p);
x2 = bpt(q);                                 
xs = linspace(x1(1),x2(1), round(abs(x1(1)-x2(1))))';
ys=linspace(x1(2),x2(2), round(abs(x1(2)-x2(2))))';
if length(xs)<=2&&length(ys)>length(xs)
    xs = interp1([x1(2) x2(2)], [x1(1) x2(1)],ys,'spline');
else
    ys = interp1([x1(1) x2(1)], [x1(2) x2(2)],xs,'spline'); 
end

curve = [curve; [xs, ys]];
if length(curve(:,1))<2
    curve=[x1;x2];
end
if ~In_nonadj_cand_pts(p,q,close_nonadj_pts,count_nonadj)
    if ~eval_conn(p,q,b,list_ind,curve)
        return
    end
end
fprintf("***Step4: Identification of Dividing Curves\n")
line=Adjust_Connect_Line(p,q,curve,Ny,Nx,b,list_ind,k,gradDX,gradDY);
line_num=line_num+1; 
Connected_Pts(line_num,1)=p;
Connected_Pts(line_num,2)=q;
conn_lines{1,line_num}=curve;


            