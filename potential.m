

clear();

scale = 6;
M = 8*scale+1.0;

N = 8*scale+1.0;
u0 = 1.0;
domX = [0 1.0];
domY = [0 1.0];




%Initialization
xGrid = linspace(domX(1),domX(2),N);
yGrid = linspace(domY(1),domY(2),M);
cells = ones(M+1,N+1);
cells(:,1)=-1.0;
cells(:,N+1)=-1.0;
cells(1,:)=0.0;
cells(M+1,:)=0.0;
    %Fill solid cells for outflow boundary
%cells(M-6*scale+1:M,(N-3*scale+1):N)=0.0;%bottom left block
%cells(M-6*scale+1:M,1:(3*scale+1))=0.0;%bottom right block
%cells(1:M/4+1,N/4+2:3*N/4+1)=0.0;%topblock  (quarters)
%cells(M/2+1:end,N/4+2:3*N/4+1)=0.0;%bottom block  (quarters)
%cells(M/2+1:M,:)=0.0;
%fill solid cells to block outflow

[solution, vx, vy] = potentialSolve(cells,domX,domY,u0);
vx = vx - ones(size(vx));
%draw cells
drawCells(domX,domY,cells);

hold on
%draw vel field
[xm, ym] = meshgrid(xGrid(2:N-1),flip(yGrid(2:M-1)));
q=quiver(xm,ym,vx,-vy);
q.LineWidth=1.0;

figure()
surf(xm,ym,log10(vx.^2+vy.^2));
colormap('jet');
colorbar();
view(2);
axis equal;
%startx= 2*dx*[1.0,1.0,1.0];
%starty= [1.0/8.0, 0.25, 3.0/8.0];

function x0 = naiiveSOR(A,b,x0,w,tol)
    nPnts=size(x0,1);
    error=1E6;
    while(error>tol)
        xp= x0;
        for i=1:nPnts
            if(i~=1)
                term=dot([A(i,1:i-1) A(i,i+1:end)],[x0(1:i-1); x0(i+1:end)]);
                %scatter(i,term)
                %hold on
            else
                term=dot(A(1,2:end),x0(2:end));
                %scatter(i,term)
                %hold on
            end
            x0(i)=(1.0-w)*x0(i)-w/A(i,i)*(b(i)-term);%sor step
            
        end
        %x0
        error=norm((x0-xp),1)
        
    end
end

function x0 = buildGuess(dx,totalPts, idLoc)
    x0 = zeros(totalPts,1);
    for i=1:totalPts
        loc= idLoc{1,i};
        x0(i)=loc(2)*dx;
    end
end

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
           row(idMap(loc(1)+1,loc(2)))=1.0;%bottom point
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
            else%corner %do backwards 1st order 2nd derivative
                row(id-2:id)=[1.0, -2.0, 1.0];
                rhs = 0.0;
             
            end
        elseif(bch==0 && bcv~=0)
            row(id-1:id+1)=row(id-1:id+1)+[1.0, -2.0, 1.0];
        end%end horizontal bc 


    %check vertical  
        if(bcv~=0)%if vertical bc
            %take care of points we know
            if(bcv>0)
%                 row(id)=row(id)-1.0;
%                 row(idMap(loc(1)-1,loc(2)))=1.0;%top point
                cellSum=sum(cells(loc(1)+2,loc(2):loc(2)+1));%sum bottom of next point
%                 curPoint = loc + [1 0];
            else
%                 row(id)=row(id)-1.0;
%                 row(idMap(loc(1)+1,loc(2)))=1.0;%bottom point
                cellSum=sum(cells(loc(1)-1,loc(2):loc(2)+1));%sum top of next point
%                 curPoint = loc + [-1 0];
            end
            if(cellSum~=0)%we are at a vertical wall
                %do backwards laplace
                row(id)=row(id)+1.0;
                row(idMap(loc(1)-bcv,loc(2)))=row(idMap(loc(1)-bcv,loc(2)))-2.0;
                row(idMap(loc(1)-2*bcv,loc(2)))=row(idMap(loc(1)-2*bcv,loc(2)))+1.0;
%                 sumLeft=sum(cells(curPoint(1):curPoint(1)+1,curPoint(2)));
%                 if(sumLeft==0)%left boundary 
%                     curPoint(2)=curPoint(2)+1;%take right point
%                 else
%                     curPoint(2)=curPoint(2)-1;%take left point
%                 end
%                 row(idMap(curPoint(1),curPoint(2)))=1.0;
%                 row(id)=row(id)-1.0;
            else
                row(idMap(loc(1)-bcv,loc(2)))=row(idMap(loc(1)-bcv,loc(2)))+1.0;
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

function drawCells(domX,domY,cells)
[M,N]=size(cells);
M=M-1;N=N-1;
dx = (domX(2)-domX(1))/(N-1);
dy = (domY(2)-domY(1))/(N-1);
cgX=linspace(domX(1)-dx,domX(2)+dx,N+2);
cgY=linspace(domY(2)+dy,domY(1)-dy,M+2);
[cx,cy]=meshgrid(cgX,cgY);
C = [[cells zeros(size(cells,1),1)] ; zeros(1,size(cells,2)+1)];
colormap('colorcube');
pcolor(cx,cy,C);
end

function [cells,idLoc,idMap,total]= createMapping(cells)
%Fills empty boundary points then ids tracked points
    [M,N]=size(cells);
    for j = 1:M
        if(cells(j,N-1)==0)%right side
            cells(j,N)=0;
        elseif(cells(j,2)==0)%left side
            cells(j,1)=0;
        end
    end
    M=M-1;
    N=N-1;
   
    idLoc = {};
    idMap = zeros(M,N);
    total=0;
    for j=2:M-1
        count=0;
        for k=2:N-1
            cellSum = sum(abs(cells(j:j+1,k))+abs(cells(j:j+1,k+1)));
            if(cellSum==4)%corner or non boundary point
                total=total+1;
                count=count+1;
                idLoc=[idLoc, [j k] ];
                idMap(j,k)=total;
            end
        end
    end
end

function [solution, vx, vy] = potentialSolve(cells,domX,domY,u0)
        %build mapping between column and grid space
    [cells,idLoc, idMap, totalPts] = createMapping(cells);
    [M,N]=size(cells);
    M=M-1;N=N-1;
    dx = (domX(2)-domX(1))/(N-1);
    dy = (domY(2)-domY(1))/(N-1);
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
    %fill array
    [A, b] = buildMatrix(cells,dx,dy,u0,uR,totalPts, idLoc, idMap);
    disp("built matrix")
    %x = A\b;
    x0 = buildGuess(dx,totalPts,idLoc);
    w=1.3;
    x = naiiveSOR(A,b,x0,w,1E-4);
    solution = rebuildGrid(x,M,N,idLoc,totalPts);
    [vx,vy] = detVel(solution,dx,M,N);
end






