

clear();


M = 12;

N = 12;
u0 = 25.0;
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

    %Fill solid cells for outflow boundary
cells(M-3:M,N-1:N)=0.0;
cells(M-3:M,1:3)=0.0;

nPnts = zeros(M,1);
nMiss = zeros(M,1);

%fill solid cells to block outflow
for j = 1:M+1
    if(cells(j,N)==0)%right side
        cells(j,N+1)=0.0;
    end
    if(cells(j,2)==0)%left side
        cells(j,1)=0.0;
    end
end


    %build mapping between column and grid space
idLoc = {};
idMap = zeros(M,N);
total=0;
for j=2:M-1
    count=0;
    for k=2:N-1
        cellSum = sum(abs(cells(j:j+1,k))+abs(cells(j:j+1,k+1)));
        if(cellSum==4 || cellSum==3 )%corner or non boundary point
            total=total+1;
            count=count+1;
            idLoc=[idLoc, [j k] ];
            idMap(j,k)=total;
        end
    end
end
totalPts = total;
rhs = zeros(totalPts);
x0 = zeros(totalPts,1);
%count missing right through cells
cLeft=0;
cRight=0;
for i=2:M
    if(cells(i,N+1)==-1)
        cRight=cRight + 1;
    end
    if(cells(i,1)==-1)
        cLeft=cLeft + 1;
    end
end

uR = u0*cLeft/cRight;%conservation motivated
%uR = u0;
%fill array
[A, b] = buildMatrix(cells,dx,dy,u0,uR,totalPts, idLoc, idMap);
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
%x0 = ones(totalPts,1);
x = sor(A,b,x0,1E-3);

solution = rebuildGrid(x,M,N,idLoc,totalPts);
[vx,vy] = detVel(solution,dx,M,N);

%draw cells
drawCells(domX,domY,dx,cells,M,N);

hold on
%draw vel field
[xm, ym] = meshgrid(xGrid(2:N-1),flip(yGrid(2:M-1)));
q=quiver(xm,ym,vx,vy);
q.LineWidth=1.0;
axis equal
streamline(xm,ym,vx,vy);



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
x=x( : ,k)
end

function [A,b] = buildMatrix(cells,dx,dy,u0,uR,totalPts, idLoc, idMap)

    A=zeros(totalPts,totalPts);
    b = zeros(totalPts,1);
    for id=1:totalPts
        loc = idLoc{1,id};
        row = zeros(1,totalPts);
        %determine if on a boundary
        bcv=0;%boundary condition flags
        bch=0;
        rhs=0.0;
        %check BC's
        %check horizontal
        shifts = [[1 0]; [0 1]; [-1 0]; [0 -1]];
        if(mod(id,totalPts/10)==0)
            disp('Filling array: (id/total)')
            disp(id/totalPts);
        end


        for i=1:4%check if there are BC's for the point
            shift=loc+shifts(i,:);
            if(idMap(shift(1),shift(2))==0)
                if(mod(i,2)==1)
                    bcv=2-i;
                else
                    bch=3-i;
                end
            end
        end

        %check for no Boundary
        if(abs(bcv)+abs(bch)==0)
           row(id-1:id+1)=[1.0,-4.0,1.0];
           row(idMap(loc(1)-1,loc(2)))=1.0;%top point
           row(idMap(loc(1),loc(2)+1))=1.0;%bottom point
           rhs=0.0;
           
        end
        %check horizontal
        if(bch~=0)%if horizontal bc
            %check if neumann or dirichlet
            if(bch>0)
                row(id-1:id)=row(id-1:id)+[1.0,-1.0];
                cellSum=sum(cells(loc(1):loc(1)+1,loc(2)+2));
            else
                row(id:id+1)=row(id:id+1)+[-1.0,1.0];
                cellSum=sum(cells(loc(1):loc(1)+1,loc(2)-1));
            end
            if(cellSum==0)%dirichlet bc
                rhs=0.0;
            elseif(cellSum<0)%through bc
                if(bch<0)
                    rhs = u0*dx;
                else
                    rhs = -uR*dx;
                end
            else%corner
                row(id)=row(id)-1.0;
                %find other referene point
                m= [0 bch];
                curPoint = loc + m;

                %decide on which way to step
                above = sum(cells(curPoint(1),curPoint(2):curPoint(2)+1));
                if(above>0 || (above==0 && sum(abs(cells(curPoint(1),curPoint(2):curPoint(2)+1)))~=0))
                    curPoint(1)=curPoint(1)-1;%move up one
                else
                    curPoint(1)=curPoint(1)+1;%move down one
                end
                %check if valid
                if(idMap(curPoint(1),curPoint(2))==0)%at a boundary point
                    %on boundary need to figure out bc of point
                    if(bch<0)
                        side=0;
                    else
                        side=1;
                    end
                    cellSum = sum(cells(curPoint(1):curPoint(1)+1,curPoint(2)+side));%sum on bc side
                    curPoint(2)=curPoint(2)-bch;%guaranteed to be this point
                    if(cellSum==0)%solid boundary
                        rhs = rhs+0.0;
                    else %through boundary
                        if(bch < 0)
                            rhs=rhs + u0*dx;
                        else
                            rhs=rhs-uR*dx;
                        end
                    end
                end   
                row(idMap(curPoint(1),curPoint(2)))=1.0;%assign value at correspoding column
            end
        elseif(bch==0 && bcv~=0)
            row(id-1:id+1)=row(id-1:id+1)+[1.0, -2.0, 1.0];
        end%end horizontal bc 


    %check vertical  
        if(bcv~=0)%if vertical bc
            %take care of points we know
            if(bcv>0)
                row(id)=row(id)-1.0;
                row(idMap(loc(1)-1,loc(2)))=1.0;%top point
                cellSum=sum(cells(loc(1)+2,loc(2):loc(2)+1));%sum bottom of next point
                curPoint = loc + [1 0];
            else
                row(id)=row(id)-1.0;
                row(idMap(loc(1)+1,loc(2)))=1.0;%bottom point
                cellSum=sum(cells(loc(1)-1,loc(2):loc(2)+1));%sum top of next point
                curPoint = loc + [-1 0];
            end
            if(cellSum~=0)%we are at a vertical wall
                sumLeft=sum(cells(curPoint(1):curPoint(1)+1,curPoint(2)));
                if(sumLeft==0)%left boundary 
                    curPoint(2)=curPoint(2)+1;%take right point
                else
                    curPoint(2)=curPoint(2)-1;%take left point
                end
                row(idMap(curPoint(1),curPoint(2)))=1.0;
                row(id)=row(id)-1.0;
            end
        elseif(bcv==0 && bch~=0)
            row(id)=row(id)-2.0;
            row(idMap(loc(1)-1,loc(2)))=1.0;%top point
            row(idMap(loc(1)+1,loc(2)))=1.0;%bottom point
        end%end veritical bc 
        A(id,:)=row;
        b(id)=rhs;
    end

end

function solution = rebuildGrid(sol,M,N,idLoc,totalPts)
  
    solution = zeros(M-2,N-2);
    
    for i=1:totalPts
        loc = idLoc{1,i};
        solution(loc(1)-1,loc(2)-1)=sol(i);
    end
            
end

function [velx,vely] = detVel(sol,dx,M,N)
    %REQUIRES CHANGE FOR GEOMETRY
    dim = size(sol);
    velx = zeros(size(sol));
    vely = zeros(size(sol));
    for i=1:dim(1)
        for j=1:dim(2)
            comy=0;
            comx=0;
            if(sol(i,j)~=0)
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
                    if(sol(i,j+1)==0)
                        velx(i,j)=(sol(i,j)-sol(i,j-1))/dx;
                    elseif(sol(i,j-1)==0)
                        velx(i,j)=(sol(i,j+1)-sol(i,j))/dx;
                    else
                        velx(i,j)=(sol(i,j+1)-sol(i,j-1))/(2.0*dx);
                    end
                end
                if(~comy)
                    if(sol(i+1,j)==0)
                        vely(i,j)=(sol(i-1,j)-sol(i,j))/dx;
                    elseif(sol(i-1,j)==0)
                        vely(i,j)=(sol(i,j)-sol(i+1,j))/dx;
                    else
                        vely(i,j)=(sol(i+1,j)-sol(i-1,j))/(2.0*dx);
                    end
                end
            end
            
        end
    end
end

function drawCells(domX,domY,dx,cells,M,N)

cgX=linspace(domX(1)-dx,domX(2)+dx,N+2);
cgY=linspace(domY(2)+dx,domY(1)-dx,M+2);
[cx,cy]=meshgrid(cgX,cgY);
C = [[cells zeros(size(cells,1),1)] ; zeros(1,size(cells,2)+1)];
colormap('colorcube');
pcolor(cx,cy,C);
end








