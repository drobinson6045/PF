global M N u0 uR totalPts nPnts nMiss;




M = 37;

N = 32;
u0 = 1.0;
domX = [0 1.0];
domY = [0 1.0];
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
nSolid = 0;%this is a bit sloppy (Has to match with filling of solid cells
    %Fill solid cells for outflow boundary
%cells(M-1:M,N-1:N)=0.0;

totalPts = (M-2)*(N-2)-nSolid;
rhs = zeros(totalPts);
A = zeros(totalPts,totalPts);
b = zeros(totalPts,1);
x0 = zeros(totalPts,1);
nPnts = zeros(M,1);
nMiss = zeros(M,1);

%fill solid cells to block outflow
for j = 1:M+1
    if(cells(j,N)==0)
        cells(j,N+1)=0.0;
    end
end


    %account for non square geometry
total=0;
for j=2:M-1
    count=0;
    for k=2:N-1
        if(sum(abs(cells(j:j+1,k))+abs(cells(j:j+1,k+1)))==4)
            total=total+1;
            count=count+1;
        end
    end
    nPnts(j) = total;
    nMiss(j)= M-2-count;
end


%count missing right cells
% count=0;
% for i=1:M
%     if(nPnts(i)-i*N == 0)
%         count=count + 1;
%     end
% end

uR = u0*(M-1)/(count-1);%conservation motivated
    

%fill array
for id=1:totalPts
    [row, rhs] = detBC(cells,id,dx,dy);
    A(id,:)=row;
    b(id)=rhs;
end
disp("built matrix")
%build initial guess of uniform flow
for i=2:M-1
   if(i==1)
       for j=1:nPnts(i)
           x0(j)=xGrid(j+1);
       end
   else
       for j=1:nPnts(i)-nPnts(i-1)
           x0(nPnts(i-1)+j)=xGrid(j+1);
       end
   end
end

x = sor(A,b,x0,1E-3);

solution = rebuildGrid(x,nPnts);
[vx,vy] = detVel(solution,dx);

%draw cells
drawCells(domX,domY,dx,cells);

hold on
%draw vel field
[xm, ym] = meshgrid(xGrid(2:N-1),yGrid(2:M-1));
quiver(xm,ym,vx,vy);




function x = sor( A ,B, x0 ,t)
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
    disp ('not convergent. ')
    disp(e)
    return
end
 

[nx, mx] = size(x0);
if nx ~= na || mx ~= 1
    disp( 'ERROR: check input')
    return
end

 
%allowed error in final answer
 

tol = t*ones(na,1);
k= 1;
 
x( : , 1 ) = x0;
err= 1E6; %intial error assumption for looping
while sum(abs(err) >= tol) ~= zeros(na,1)
    x( : ,k+ 1 ) = (D+w*L) \ ( D*(1-w) - w*U)*x( : ,k) + (D+ w*L)\B;% SOR formula
    err = x( :,k+1) - x( :, k);% finding error
    k = k + 1; 
end
 
fprintf (' The final ans obtained after %d iterations\n', k)
x( : ,k)
end

function [row, rhs] = detBC(cells,id,dx,dy)
    global M N u0 uR totalPts nPnts nMiss;
    rr = find((nPnts-id)>=0,1);
    if(rr>1)
        cc = 1+id-nPnts(rr-1);
    end
    row = zeros(1,totalPts);
    %determine if on a boundary
    bcv=0;%boundary condition flags
    bch=0;
    if(cc==2)
        bch=-1;
    elseif(cc==N-1)
        bch=1;
    end
    if(rr==2)
        bcv=1;
    elseif(rr==M-1)
        bcv=-1;
    end
    %check BC's
    %check for no Boundary
    if(abs(bch)+abs(bcv)==0)
       row(id-1:id+1)=[1.0,-4.0,1.0];
       row(id-(M-2-nMiss(rr-1)))=1.0;%top point
       row(nPnts(rr)+cc-1)=1.0;%bottom point
       rhs=0.0;
       return
    end
    %check horizontal
    if(bch==-1)%left side through
        row(id:id+1)=row(id:id+1)+[-1.0,1.0];
        rhs = u0*dx;
    elseif(bch==1)%right side through
        row(id-1:id)=row(id-1:id)+[1.0,-1.0];
        rhs = -u0*dx;
    elseif(bch==0)
        row(id-1:id+1)=[1.0,-2.0,1.0];
        rhs=0.0;
    end
    %check vertical
    if(bcv==-1)%solid bottom
        row(id)=row(id)-1.0;
        row(id-(M-2-nMiss(rr-1)))=row(id-(M-2-nMiss(rr-1)))+1.0;
    elseif(bcv==1)%top solid
        row(id) = row(id)-1.0;
        row(nPnts(rr)+cc-1)=row(nPnts(rr)+cc-1)+1.0;
    elseif(bcv==0)
        row(id-(M-2-nMiss(rr-1)))=1.0;%top point
        row(nPnts(rr)+cc-1)=1.0;%bottom point
        row(id)=row(id)-2.0;
    end

end

function solution = rebuildGrid(sol,nPnts)
    global M N;
    solution = zeros(M-2,N-2);
    id = 1;
    for i=2:M-1
        for j =1:(nPnts(i)-nPnts(i-1))
            solution(i-1,j)=sol(id);
            id=id+1;
        end
    end
            
end

function [velx,vely] = detVel(sol,dx)
    %REQUIRES CHANGE FOR GEOMETRY
    global M N;
    dim = size(sol);
    velx = zeros(size(sol));
    vely = zeros(size(sol));
    for i=1:dim(1)
        for j=1:dim(2)
            comy=0;
            comx=0;
            if(i==1)%top boundary
                vely(i,j)=(sol(i,j)-sol(i+1,j))/dx;
                comy=comy+1;
            end
            if(i==M-2)
                vely(i,j)=(sol(i-1,j)-sol(i,j))/dx;
                comy=comy+1;
            end
            if(j==1)
                velx(i,j)=(sol(i,j+1)-sol(i,j))/dx;
                comx=comx+1;
            end
            if(j==N-2)
                velx(i,j)=(sol(i,j)-sol(i,j-1))/dx;
                comx=comx+1;
            end
            if(~comx)
                velx(i,j)=(sol(i,j+1)-sol(i,j-1))/(2.0*dx);
            end
            if(~comy)
                vely(i,j)=(sol(i+1,j)-sol(i-1,j))/(2.0*dx);
            end
            
        end
    end
end

function drawCells(domX,domY,dx,cells)
global N M;
cgX=linspace(domX(1)-dx,domX(2)+dx,N+2);
cgY=linspace(domY(1)-dx,domY(2)+dx,M+2);
[cx,cy]=meshgrid(cgX,cgY);
C = [[cells zeros(size(cells,1),1)] ; zeros(1,size(cells,2)+1)];
colormap('colorcube');
pcolor(cx,cy,C);
end








