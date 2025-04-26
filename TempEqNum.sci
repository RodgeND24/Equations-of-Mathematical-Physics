clc; 
//функция численного решения уравнения теплопроводнойсти
function [u,x,t] = TempEqNum(Nx,Nt,L,T,a)
h=L/Nx;
tau=T/Nt;


for i=1:Nx+1
    x(i)=(i-1)*h;
    u(i,1)=u_0(x(i));
end

for j=1:Nt+1
    t(j)=(j-1)*tau;
    u(1,j)=u_1(t(j));
    u(Nx+1,j)=u_2(t(j));
end

gam=a^2*tau/h^2;

for j=1:Nt
    for i=2:Nx
        //u(i,j+1)=gam*(u(i+1,j)-2*u(i,j)+u(i-1,j))+2*u(i,j)-u(i,j-1);
        u(i,j+1)=gam*u(i+1,j)+(1-2*gam)*u(i,j)+gam*u(i-1,j)+tau*f(x(i),t(j));
    end
end
endfunction

//Функция в правой части уравнения
function func = f(x,t)
    func = 0
endfunction

//Граничные условия
function u1 = u_1(t)
    u1 = 0
endfunction

function u2 = u_2(t)
    u2 = 0
endfunction

//начальные условия
function u0 = u_0(x)
    //u0 = 2*exp(-(x-2)^2)
    u0 = 100*sin(%pi*x/2)^2
endfunction


L=4;
T=60;
a=0.1;
Nx=40;
Nt=2*(a*Nx/L)^2*T;
disp(T/Nt);
[U,X,T]=TempEqNum(Nx,Nt,L,T,a);

surf(X,T,U');
xlabel('X');
ylabel('T');
zlabel('U');

// Сохранение результата в формате .gif

Max=max(max(U));
Min=min(min(U));
h_fig = figure; 
h_fig.background = 8;
h_point = plot(X, U(:,1));
xgrid();
gca().data_bounds = [0, L, Min, Max];


outgif = 'Temp.gif';
mdelete(outgif);
idGif = animaGIF(gcf(), outgif, 10, 0);

for j=2:(length(T+1)*1)
    drawlater(); 
    p=int((j-1)./1)+1;
    h_point.data = [X,U(:,p)];
    drawnow();
    idGif = animaGIF(gcf(), idGif);
end
animaGIF(idGif);
clear;


