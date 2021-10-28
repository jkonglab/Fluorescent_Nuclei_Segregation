function Nuclei_Seg(SIGMA,ARC,QUALTHRESH,QUAL,NumberOfimages,Min_Distance,Max_Distance,ANGLE,UNUSED,Val_Nonclose)

% Nuclei_Seg performs clumped nuclei segmentation
% Params:
%   SIGMA is the value used to define gaussian kernel, to obtain the
%   image derivates
%   ARC is the arc length threshold used to identify whether two adjacent
%   points are close or not
%   QUALTHRESH is the ellipse fitting threshold
%   QUAL is also an ellipse fitting threshold which is used when there is
%   no good candidates selected, QUAL < QUALTHRESH 
%   Numberofimages is the number of image patches used to obtain the
%   average nucleus size
%   Min_Distance is r_1 defined in paper "Clumped Nuclei Segmentation with
%   Adjacent Point Match and Local Shape-Based Intensity Analysis in Fluor-
%   escence Microscopy Images"
%   Max_Distance is r_2 defined in the same paper 
%   ANGLE is the angle threshold of two segmentation lines
%   UNUSED is used to decide whether there are enough points that should be
%   splitted futher
%   Vak_Nonclose is the value threshold of two nonadjcent points 
%   =======================================================================================
%   Copyright (C) 2018  Xiaoyuan Guo
%   Email: xiaoyuan.guo@emory.edu
%   =======================================================================================

if nargin==0 
    SIGMA = 6;
    ARC = 30;
    QUALTHRESH= 0.70;
    QUAL= 0.60;
    NumberOfimages=10;
    Min_Distance=45;%30
    Max_Distance=76;%65
    ANGLE=80;
    UNUSED=1;
    Val_Nonclose=200;
    Nuclei_Seg(SIGMA,ARC,QUALTHRESH,QUAL,NumberOfimages,Min_Distance,Max_Distance,ANGLE,UNUSED,Val_Nonclose);
end

%% global parameters
warning off;
global Min_Distance;
global Max_Distance;
global QUALTHRESH;
global QUAL;
global blue;
global avg_area;
global bw_poly;
global gradDX1;
global gradDY1;

avg_area=Avg_Area_Nucleus(NumberOfimages);
%% image pre-processing
% I=imread('/Users/xiaoyuanguo/Desktop/Data/test_nuclei7.tif');
I=imread('/Users/xiaoyuanguo/Desktop/FISH/1_S16-28151_A2/data/S16_28151 63x_BT_3.tif');
blue = I(:,:,3);
bw=im2bw(blue,0.25);
bw_cleaned=bwareaopen(bw,60);
bw_cleaned = imfill(bw,'holes');
[~, L] = bwboundaries(bw_cleaned);
stats = regionprops(L, 'Area');
bw_cleaned = bwareaopen(bw_cleaned, round(mean([stats.Area])/10));
[B, L] = bwboundaries(bw_cleaned);
stats = regionprops(L, 'Area','Centroid');

%% find high curvature points

% define gaussian kernel and its derivatives
t = ceil(3*SIGMA);
h = exp(-(-t:t).^2/(2*SIGMA^2)) / (SIGMA*sqrt(2*pi));
h = h/sum(h);

h1 = conv(h, [1 0 -1], 'valid');
h1 = h1 - sum(h1)/length(h1); %VERY IMPORTANT TO BRING sum(h1) to ZERO
h2 = conv(h1, [1 0 -1], 'valid');
h2 = h2 - sum(h2)/length(h2); %VERY IMPORTANT TO BRING sum(h2) to ZERO

imshow(I);
I = mat2gray(blue);
I = imgaussfilt(I, 2);

[Gmag, Gdir] = imgradient(I,'prewitt');
hold on;

%% Segmentation clumped nuclei one by one

for i = 1:length(B)
    fprintf("Segment clump nuclei %d\n",i);
    fprintf("***************************\n\n")
    b = B{i};
    ar=area(b);
    
    %plot clumped nuclei boundaries 
    boundary=B{i};
    plot(boundary(:,2),boundary(:,1),'r','linewidth',3);
   
    b(:,1) = imfilter(b(:,1), h' , 'circular', 'same', 'conv');
    b(:,2) = imfilter(b(:,2), h' , 'circular', 'same', 'conv');
    plot(b(:,2), b(:,1), 'w', 'LineWidth', 2);
    
    % convolution with 1st and 2nd derivative of gaussian filter
    x1 = imfilter(B{i}(:,2), h1', 'circular', 'same', 'conv'); % conv2(B{i}(:,2), h1,'same') also works
    y1 = imfilter(B{i}(:,1), h1', 'circular', 'same', 'conv');
    x2 = imfilter(B{i}(:,2), h2', 'circular', 'same', 'conv');
    y2 = imfilter(B{i}(:,1), h2', 'circular', 'same', 'conv');
    %eps for return the minimum point precision
    %compute curvature values of every boundary point
    k = (x1.*y2 - y1.*x2)./(eps+(x1.^2+y1.^2).^(3/2));
    
    %% High curvature candidate point voting
    fprintf("***Step1: High curvature points voting\n")
    [curvature,list_ind,n]=Find_Curv_Candidate_Points(k,b); 
    
    arclens = arrayfun(@(x1,x2,y1,y2) sqrt((x1-x2)^2+(y1-y2)^2), b(:,1), b([end 1:end-1],1),b(:,2), b([end 1:end-1],2));

    if n>1
    %label the candidate points 
        plot(b(list_ind,2), b(list_ind,1), 'yo', 'MarkerSize', 10);
        for t = 1:length(list_ind)
            str=sprintf('%2d',t);
            text(b(list_ind(t),2),b(list_ind(t),1),str,'fontsize',20,'color',[1,1,1]);
        end  
        
        % When the area of the cluster is too small, return
        if  ar<1.5*avg_area
            continue;   
        end
        
        candidate_pts=linspace(1,n,n);%generate points index
        fprintf("***Step2:Detect close Adjacent points and Nonadjacent points using radius r_1\n")
        %% compute Euclidean distance between two points
        close_Mat=Close_Detection(list_ind,b);%close_Mat stores the distance value of each two points 
        
        %% Detect close pointpairs, close_adj_pts stores the adjacent
        %pointpairs, close_nonadj_pts stores nonadjacent pointpairs,
        %count_nonadj records the number of close nonadjacent pointpairs
        [close_adj_pts,close_nonadj_pts,count_nonadj,count_adj,nonadj_close_distance]=close_pts(close_Mat,curvature,list_ind,b,Gmag, Gdir);
        
        Connected_Pts=zeros(2*n,2);%store connection information
        line_num=0;% record the number of connected lines
        conn_lines={1,n};% Init the conn_lines set
        %We split the whole boundary with close nonadjacent pointpairs 
        Candidate_pts={};%store the candidate pts for each seperated boundary
        % to filter those point pairs that don't need to merge
        particles = b(list_ind,:);
        Np = size(particles,1);
        bounds = cell(Np,1);
        [Ny, Nx, ~] = size(I);
        %convert region of interest polygon to region mask
        bw_poly = poly2mask(b(:,2), b(:,1), size(I,1),size(I,2));
        bw_poly = imdilate(bw_poly, ones(5,5));    
        D = graydist(double(I), ~bw_poly, 'quasi-euclidean');
        [gradDX,gradDY] = gradient(double(D),0.2,0.2);
        L2 = sqrt(gradDX.^2 + gradDY.^2);
        gradDX1=gradDX;
        gradDY1=gradDY;
        gradDX = gradDX./L2;
        gradDY = gradDY./L2;
        fprintf("***Step2: Close Point Pair Screening\n");
        fprintf("*******Step2.1: Analysis for Pairs of Adjacent Points in Proximit\n")
        [close_adj_pts,count_adj,conn_lines,Connected_Pts,line_num]=Augment_adj_Close(close_adj_pts,count_adj,list_ind,k,b,Gdir,conn_lines,Connected_Pts,line_num,Nx,Ny,gradDX,gradDY,close_nonadj_pts,count_nonadj);     
        remove_list=[];
        fprintf("*******Step2.2: Analysis for Pairs of Non-Adjacent Points in Proximity\n")
        if count_adj>0
            [candidate_pts,curvature,new_close_adj_pts,new_count_adj,close_nonadj_pts,count_nonadj]=Merge_Close_Adj_Pts(b,close_adj_pts,curvature,list_ind,close_nonadj_pts,count_nonadj,count_adj,nonadj_close_distance,candidate_pts,remove_list);            
        elseif count_nonadj>0
            [close_nonadj_pts,count_nonadj]=Update_Close_Nonadj_pts(close_adj_pts,count_adj,close_nonadj_pts,count_nonadj,remove_list,nonadj_close_distance);
           
        end
%         bx = b(:,2);
%         by = b(:,1);
       
        % divergence
        % [gX, gY] = imgradientxy(255-I,'central');
        % div = divergence(gX,gY);

        %% dynamic gradients
%         xidx = 1:5:Nx;
%         yidx = 1:5:Ny;
% 
%         [X,Y] = ndgrid(xidx,yidx);
%         [dx dy]=gradient(I)    
%        
%          plot(bx, by, 'w', 'LineWidth', 2);
%          qu=quiver(X', Y',gradDX(yidx,xidx),gradDY(yidx,xidx));
%          c = qu.Color;
%          qu.Color = 'green';
%          hold on
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Split Candidate Points(Begin)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_Candi_sets=1;
        if count_nonadj>0
            Candidate_pts={1,length(candidate_pts)+1};
            for j=1:count_nonadj
                p=close_nonadj_pts(j,1);
                q=close_nonadj_pts(j,2);
                if j==1||(j>1&&line_num==0)
                    last_line_num=line_num;
                    [conn_lines,Connected_Pts,line_num]=Connect_pts(p,q,b,list_ind,conn_lines,Connected_Pts,line_num,Nx,Ny,k,gradDX,gradDY,close_nonadj_pts,count_nonadj);
                    if last_line_num==line_num
                        continue;
                    end
                    [cand_pts1,cand_pts2]=Split_Candidate_Pts(p,q,candidate_pts);
                    Candidate_pts{1,1}=cand_pts1;
                    Candidate_pts{1,2}=cand_pts2;
                    num_Candi_sets=num_Candi_sets+1;
                else
                    [p,q]=Swap_pq(p,q);
                    is_conn=line_num;
                    [conn_lines,Connected_Pts,line_num]=Connect_pts(p,q,b,list_ind,conn_lines,Connected_Pts,line_num,Nx,Ny,k,gradDX,gradDY,close_nonadj_pts,count_nonadj);
                    is_conn=abs(is_conn-line_num);
                    if is_conn==0
                        continue
                    end

                    [new_candidate_pts,where]=In_Where_Candidatet(p,q, Candidate_pts,num_Candi_sets);
                    if ~where&&line_num
                          continue
                    end
                     
                    [cand_pts1,cand_pts2]=Split_Candidate_Pts(p,q,new_candidate_pts);
                    num_Candi_sets=num_Candi_sets+1;
                    Candidate_pts{1,where}=cand_pts1;
                    Candidate_pts{1,num_Candi_sets}=cand_pts2;                       
                end
            end
        
        end
%%%%%%%%%%%%%%Split Candidate Points(End)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf("***Step 3: Point Pair Segmentation Assessment by Ellipse Fitting\n");
        for it=1:num_Candi_sets
            if num_Candi_sets==1
                cann_pts=candidate_pts;           
            else 
            %evaluate the situation, if point number is small and the area
            %is small, stop segmentation
                ret=Eval_Cann_pts_Seg(Candidate_pts{1,it},b,list_ind,close_nonadj_pts,Connected_Pts,conn_lines,count_nonadj);
                if length(Candidate_pts{1,it})<=3
                    continue
                elseif ret
                    continue
                elseif length(Candidate_pts{1,it})==4
                    cann_pts_temp=Candidate_pts{1,it};
                    [ret1,num]=In_Connected_Pts(cann_pts_temp(1),cann_pts_temp(2),Connected_Pts);
                    [ret2,num]=In_Connected_Pts(cann_pts_temp(3),cann_pts_temp(4),Connected_Pts);
                    [ret3,num]=In_Connected_Pts(cann_pts_temp(1),cann_pts_temp(4),Connected_Pts);
                    [ret4,num]=In_Connected_Pts(cann_pts_temp(2),cann_pts_temp(3),Connected_Pts);
                if (ret1&&ret2)||(ret3&&ret4)
                    continue
                end 
            end           
            cann_pts= Candidate_pts{1,it};      
            end      
       [good_conn_pts_p,good_conn_pts_q,count]=Compute_Good_adj_Qual(list_ind,b,cann_pts,conn_lines,Connected_Pts);       
       go_on=true;

       if count 
           unused_num=length(cann_pts);
            for j=1:count%count is the number of good pts
                if (unused_num<UNUSED&& go_on==false)
                     continue
                end
                new_p=good_conn_pts_p(j);
                new_q=good_conn_pts_q(j);
                last_line_num=line_num;
                [conn_lines,Connected_Pts,line_num]=Connect_pts(new_p,new_q,b,list_ind,conn_lines,Connected_Pts,line_num,Nx,Ny,k,gradDX,gradDY,close_nonadj_pts,count_nonadj);
                if last_line_num==line_num
                    continue
                end
                [unused_num,~,go_on]=Left_Unconnected(cann_pts,Connected_Pts);
                 
            end
       end
        
        %use count we could find the last filled connect lies
        [New_good_pts_num, New_good_Ps,New_good_Qs]=Find_New_Good_Pts(list_ind,b,cann_pts,conn_lines,line_num,Connected_Pts,count,close_nonadj_pts,count_nonadj);
        %To control the same point ptpairs repeated
        
        while New_good_pts_num       
            for t=1:New_good_pts_num
                st=false;
                new_p=New_good_Ps(t);
                new_q=New_good_Qs(t);
                last_line_num=line_num;
                [conn_lines,Connected_Pts,line_num]=Connect_pts(new_p,new_q,b,list_ind,conn_lines,Connected_Pts,line_num,Nx,Ny,k,gradDX,gradDY,close_nonadj_pts,count_nonadj); 
                [unused_num,~,go_on]=Left_Unconnected(cann_pts,Connected_Pts);
                 if unused_num<=UNUSED|| ~go_on || (line_num==last_line_num)
                     st=true;
                    break;
                 end 
            end
             if unused_num<=UNUSED|| ~go_on ||st
                 break;
             end            
                [New_good_pts_num, New_good_Ps,New_good_Qs]=Find_New_Good_Pts(list_ind,b,cann_pts,conn_lines,line_num,Connected_Pts,count,close_nonadj_pts,count_nonadj);           
        end          
    end
  end
end

hold off;


