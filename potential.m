global M N u0 uR totalPts nPnts nMiss;




M = 50;

N = 50;
u0 = 2.0;
domX = [0 100.0];
domY = [0 100.0];
dx = (domX(2)-domX(1))/(N-1);
dy = (domY(2)-domY(1))/(N-1);
xGrid = linspace(domX(1),domX(2),N);
yGrid = linspace(domY(1),domY(2),M);
cells = ones(M+1,N+1);
cells(:,1)=-1.0;
cells(:,N+1)=-1.0;
cells(1,:)=0.0;
cells(M+1,:)=0.0;


%Initialization
nSolid = 4;%this is a bit sloppy (Has to match with filling of solid cells
    %Fill solid cells for outflow boundary
cells(M-1:M,N-1:N)=0.0;

totalPts = M*N-nSolid;
rhs = zeros(totalPts);
A = zeros(totalPts,totalPts);
b = zeros(totalPts,1);
x0 = zeros(totalPts,1);
nPnts = M*ones(M,1);
nMiss = zeros(M,1);
%fill solid cells
for j = 1:M+1
    if(cells(j,N)==0)
        cells(j,N+1)=0.0;
    end
end


    %account for non square geometry
total=0;
for j=1:M
    count=0;
    for k=1:N
        if(sum(abs(cells(j:j+1,k))+abs(cells(j:j+1,k+1)))~=0)
            total=total+1;
            count=count+1;
        end
    end
    nPnts(j) = total;
    nMiss(j)= M-count;
end


%count missing right cells
count=0;
for i=1:M
    if(nPnts(i)-i*N == 0)
        count=count + 1;
    end
end

uR = u0*(M-1)/(count-1);%conservation motivated
    

%fill array
for id=1:totalPts
    [row, rhs] = detBC(cells,id,dx,dy);
    A(id,:)=row;
    b(id)=rhs;
end
disp("built matrix")
%build initial guess of uniform flow
for i=1:M
   if(i==1)
       for j=1:nPnts(i)
           x0(j)=xGrid(j);
       end
   else
       for j=1:nPnts(i)-nPnts(i-1)
           x0(nPnts(i-1)+j)=xGrid(j);
       end
   end
end

x = SOR(A,b,x0,1E-3);

function X = SOR( A ,B, X0 ,t)
[na , ma ] = size (A);
[nb , mb ] = size (B);
if nb ~= na || mb~=1
   disp( 'ERROR:input error..pls re-enter data')
   return
end
 
w=1.3;
% A = D + L + U
 
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
 
% check for convergence
e= max(eig(( D+w*L) \ ( D*(1-w) - w*U)))
 
if abs(e) >= 1
    disp ('largest eigen value of iterative matrix is not less than 1') 
    disp ('not convergent. Please try some other process.')
    disp(e)
    return
end
 

[nx, mx] = size(X0);
if nx ~= na || mx ~= 1
    disp( 'ERROR: pls check input')
    return
end

 
%allowed error in final answer
 

tol = t*ones(na,1);
k= 1;
 
X( : , 1 ) = X0
err= 1E6; %intial error assumtion for looping
while sum(abs(err) >= tol) ~= zeros(na,1)
    X ( : ,k+ 1 ) = (D+w*L) \ ( D*(1-w) - w*U)*X( : ,k) + inv(D+ w*L)*B;% SOR formula
    err = X( :,k+1) - X( :, k);% finding error
    disp(err)
    k = k + 1; 
end
 
fprintf (' The final ans obtaibed after %d itaration which is  \n', k)
X( : ,k)
end








function [row, rhs] = detBC(cells,id,dx,dy)
    global M N u0 uR totalPts nPnts nMiss;
    rr = find((nPnts-id)>=0,1);
    if(rr>1)
        cc = id-nPnts(rr-1);
    else
        cc = id;
    end
    row = zeros(1,totalPts);
    grid = cells(rr:rr+1,cc:cc+1);
    
    if(sum(grid(1,:)+grid(2,:))==4 ||sum(grid(1,:)+grid(2,:))==3)%if an interior point no BCS just laplace
        row(id-1:id+1)=[1.0, -4.0, 1.0]/(dx*dx);
        if(rr>=2)
            row(id-(M-nMiss(rr-1)))=1.0/(dy*dy);%top point
            row(nPnts(rr)+cc)=1.0/(dy*dy);%bottom point
        %else
        %    row(cc)=1.0;
        %    row(nPnts(rr)+cc)=1.0;
        end
        rhs=0.0;
    else
        if(sum(abs(grid(:,1)))==0)%left solid bc
            %none yet algorithm can't catch this yet
        end
        if(sum(abs(grid(:,2)))==0)%right solid bc
            row(id-1:id)=row(id-1:id)+[-1.0 1.0]/dx;
            row(id-2:id)=row(id-2:id)+[1.0 -2.0 1.0]/(dx*dx);  %backward laplace
            
        end
        if(sum(abs(grid(1,:)))==0)%top solid bc
            row(id) = row(id)+1.0/dy;
            row(nPnts(rr)+cc)=row(nPnts(rr)+cc)-1.0/dy;%these might need sign change (depend on nhat)
            %backward Laplace (y)
            row(id)=row(id)+1.0/(dy*dy);
            row(nPnts(rr)+cc)=row(nPnts(rr)+cc)-2.0/(dy*dy);
            row(nPnts(rr+1)+cc)=row(nPnts(rr)+cc)+1.0/(dy*dy);      
        end
        if(sum(abs(grid(2,:)))==0)%bottom solid bc
            row(id)=row(id)-1.0/dy;
            row(id-(M-nMiss(rr-1)))=row(id-(M-nMiss(rr-1)))+1.0/dy;
            %forward Laplace (y)
            row(id)=row(id)+1.0/(dy*dy);%center point
            row(id-(M-nMiss(rr-1)))=row(id-(M-nMiss(rr-1)))-2.0/(dy*dy);
            row(id-(2*M-nMiss(rr-1)-nMiss(rr-2)))=row(id-(2*M-nMiss(rr-1)-nMiss(rr-2)))+1.0/(dy*dy);
        end
        if(sum(grid(:,2))<=0)%right through bc
           row(id-1:id)=row(id-1:id)+[-1.0 1.0]/dx; 
           row(id-2:id)=row(id-2:id)+[1.0 -2.0 1.0]/(dx*dx);  %backward laplace
           rhs=uR;
        end
        if(sum(grid(:,1))<=0)%left through bc
           row(id:id+1)=row(id:id+1)-[-1.0 1.0]/dx;%might need sign change
           row(id:id+2)=row(id:id+2)+[1.0 -2.0 1.0]/(dx*dx);  %backward laplace
           rhs=u0;%rhs needed
        end  
        %%Add forward backward or centered laplace to cells x or y laplace
        if((sum(grid(:,1))==1 || sum(grid(:,1))==2)&& (sum(grid(:,2))==1 || sum(grid(:,2))==2))%no left or right BC
            row(id-1:id+1)=row(id-1:id+1)+[1.0, -2.0, 1.0]/(dx*dx);
            rhs=0.0;
        end
        if((sum(grid(1,:))==1 || sum(grid(1,:))==2)&& (sum(grid(2,:))==1 || sum(grid(2,:))==2))%no top or bottom BC (corner)
            row(id)=row(id)-2.0/(dy*dy);%center point
            row(id-(M-nMiss(rr-1)))=row(id-(M-nMiss(rr-1)))+1.0/(dy*dy);
            row(nPnts(rr)+cc)=row(nPnts(rr)+cc)-1.0/(dy*dy);
            rhs=0.0;
        end
        
    end
end









