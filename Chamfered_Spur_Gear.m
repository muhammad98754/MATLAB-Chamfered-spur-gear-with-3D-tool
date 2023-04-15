%Drawing Chamfered gear
m_n=5
z=50
alpha_n=20*pi/180
x_coef=0
h_aP_coef=1+x_coef
h_fP_coef=1.25+x_coef
k=0
b=5
r_p=(m_n*z)/2
r_b=r_p*cos(alpha_n)
r_a= z*m_n/2 + m_n*(h_aP_coef + x_coef + k);
r_f= z*m_n/2 - m_n*(h_fP_coef - x_coef)
r_y=r_b:0.1:r_a
alpha_y=acos(r_b./r_y)
%gear profile involute curve  
x=r_y.*sin(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_y)+alpha_y))+(2*x_coef*tan(alpha_n)*m_n)
y=r_y.*cos(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_y)+alpha_y))+(2*x_coef*tan(alpha_n)*m_n)
plot(x,y)
hold on
%opposite side of gear profile
mirror=[flip(-x); flip(y)]
plot(mirror(1,:),mirror(2,:))
% Top side gear profile
x_a=x(1,end);
y_a=y(1,end);
theta=atan(x_a./y_a);
r=sqrt((x(1,end))^2+(y(1,end))^2)
t=theta:-0.01:-theta
x_t=r.*sin(t)
y_t=r.*cos(t)
plot(x_t,y_t)
% bottom land gear profile and application of fillet radius
dy = cos(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_y)+alpha_y))
dx = sin(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_y)+alpha_y))
r_c=abs(r_f-r_b)
k = dy/dx
theta_f = atan(-(1/k))
x_start = x(1,1)
y_start = y(1,1)
X_c = x_start + (r_c*cos(theta_f))
Y_c = y_start + (r_c*sin(theta_f))
x_f = linspace(x_start,X_c,150)
%eqn of a circle: (x-a)^2 + (y-b)^2 = r^2
y_f = Y_c - (sqrt(r_c^2 - (x_f-X_c).^2))
 y_f(1,1) = y_start
plot(x_f,y_f,"G")
plot(flip(-x_f),flip(y_f),"G")
a=[flip(x_f),x,mirror(1,:),-x_f]
s=[flip(y_f),y,mirror(2,:),y_f]
[x_tooth, y_tooth]=xy(a,s)
% chamfer parametrs
delta_a= 35*pi/180;
a_ch= +1;
d_ch_flat= 2*r_f-3 ;
d_ch_down= 2*r_f-1 ;
r_e= 1.2;
% convert diameter into radius
r_ch_flat= d_ch_flat/2;
r_ch_down= d_ch_down/2;
% Insert the point for d_ch_down
n_ch_down= find(r_ch_down < sqrt(x_tooth.^2 + y_tooth.^2),1,'first');
[x_tooth, y_tooth]= xy_ch_down(x_tooth, y_tooth, r_ch_down, n_ch_down);
% Add point end to the chamfer at tip diameter
n_ch_a= find(x_tooth < 0,1,'first');
x_tooth= [x_tooth(1:n_ch_a-1),a_ch,x_tooth(n_ch_a:end)];
y_tooth= [y_tooth(1:n_ch_a-1),sqrt(r_a^2 - a_ch^2),y_tooth(n_ch_a:end)];
% Convert profile into vertices two profiles on 0 and -b height
V_tooth= [x_tooth',y_tooth', zeros(size(x_tooth')); ...
          x_tooth',y_tooth', -b.*ones(size(x_tooth')) ];
% Change vertices in flat chamfer area     
V_tooth(n_ch_down:n_ch_a-1,3)= -tan(delta_a)*(V_tooth(n_ch_down:n_ch_a-1,1)-a_ch);
% Make faces for periphery
ns_top= 1:size(x_tooth');
ns_bottom= size(x_tooth')+1:size(V_tooth,1);
F_periphery= [ns_bottom(1:end-1)',ns_bottom(2:end)',ns_top(2:end)',ns_top(1:end-1)'];
% Add point ch_flat to vertices
p_ch_flat= [a_ch, sqrt(r_ch_flat^2 - a_ch^2), 0];
V_tooth= [V_tooth;p_ch_flat];
n_ch_flat= size(V_tooth,1);
% make flat chamfer face
F_ch_flat= [n_ch_down:n_ch_a,n_ch_flat];
p_ch_down= V_tooth(n_ch_down,:);
% Move and rotate the vertices so that the line from p_ch_flat to p_ch_down
% is in line with the y-axis!
V_tooth= V_tooth - p_ch_flat;
V_tooth= Mult_homoMat(makehgtform('yrotate',-delta_a), V_tooth);
phi_z= atan2(V_tooth(n_ch_down,1),V_tooth(n_ch_down,2));
V_tooth= Mult_homoMat(makehgtform('zrotate', phi_z), V_tooth);
V_tooth(:,3)= V_tooth(:,3)- r_e;
% rotate the vector [0,0,1] z-axis with the roation on the vertices
vec001= Mult_homoMat(makehgtform('zrotate', phi_z)* ...
                     makehgtform('yrotate',-delta_a), [0,0,1]);
% calulation of the intersection between cylinder by r_e and perphery
ns_corner_down=[];
for n=n_ch_down:-1:1
    [p_out,intersec]= Cylinder2line(r_e, V_tooth(n,:), vec001);
    if intersec==true
       V_tooth(n,:)= p_out;
       ns_corner_down= [ns_corner_down,n];
    end
end
% Move and rotate the vertices back into initial state
V_tooth(:,3)= V_tooth(:,3)+ r_e;
V_tooth= Mult_homoMat(makehgtform('zrotate', -phi_z), V_tooth);
V_tooth= Mult_homoMat(makehgtform('yrotate', delta_a), V_tooth);
V_tooth= V_tooth + p_ch_flat;
% Calculate the intesection of the cylinder r_e with x/y plane
v_ch_flat2ch_down= p_ch_flat-p_ch_down;
ns_corner_top= [1:length(ns_corner_down)]+ size(V_tooth,1);
V_tooth= [V_tooth;-V_tooth(ns_corner_down,3)/v_ch_flat2ch_down(3).*ones(size(ns_corner_down))'*v_ch_flat2ch_down + V_tooth(ns_corner_down,:)]
% make the faces
F_corner= [ns_corner_down(1:end-1)',ns_corner_top(1:end-1)',ns_corner_top(2:end)',ns_corner_down(2:end)'];
F_triangle=[ns_corner_down(end),ns_corner_top(end),ns_corner_down(end)-1];
F_2D=[1:ns_corner_down(end)-1,flip(ns_corner_top),n_ch_a:ns_top(1,end)]
gear=[]
clf
axis equal
for i=0:2*pi/z:2*pi
    V_tooth_temp=makehgtform("zrotate",i)*[V_tooth';ones(1,length(V_tooth))]
    V_tooth_rotated=[V_tooth_temp(1,:)', V_tooth_temp(2,:)' V_tooth_temp(3,:)']
    V_tooth_b=makehgtform("zrotate",i)*[a;s;-b.*ones(1,length(a));ones(1,length(a))]
    g1=patch('Vertices',V_tooth_rotated,'Faces',F_ch_flat,'FaceColor','Yellow')
    g2=patch('Vertices',V_tooth_rotated,'Faces',F_corner,'FaceColor','Yellow')
    g3=patch('Vertices',V_tooth_rotated,'Faces',F_triangle,'FaceColor','yellow')
    g4=patch('Vertices',V_tooth_rotated,'Faces',F_periphery,'FaceColor',[0.7,0.7,0.7])
    g5=patch('Vertices',V_tooth_rotated,'Faces',F_2D,'FaceColor','blue')
    g6=patch(V_tooth_b(1,:),V_tooth_b(2,:),V_tooth_b(3,:),'Yellow')
    hold on
    plot3(V_tooth(:,1),V_tooth(:,2),V_tooth(:,3))
    concat_gear=[g1;g2;g3;g4;g5;g6]
    gear=[gear;concat_gear]
end
hold on

r_root=sqrt((x_f(1,end))^2+(y_f(1,end))^2)
t_root=0:0.01:2*pi
x_b=r_root.*sin(t_root)
y_b=r_root.*cos(t_root)
plot(x_b,y_b,'green')
z_b=horzcat(zeros(1,length(x_b)),-b.*ones(1,length(x_b)))
g7=patch(x_b,y_b,zeros(1,length(x_b)),'green')
g8=patch(x_b,y_b,-b.*ones(1,length(x_b)),'green')
g9=patch([x_b,x_b],[flip(y_b),y_b],z_b,'green')
view(3)


xlim([-200 200])
ylim([-200 200])
zlim([-10 10])

function [x_tooth, y_tooth]= xy_ch_down(x_tooth, y_tooth, r_ch_down, n_ch_down)
% Calculates with a 4 point cubic spline the point p_ch_down
% add the point p_ch_down to the profile
   [phi_tooth,r_tooth]= cart2pol(x_tooth(n_ch_down-2:n_ch_down+1),y_tooth(n_ch_down-2:n_ch_down+1));
   phi_ch_down=spline(r_tooth, phi_tooth, r_ch_down);
   x_tooth= [x_tooth(1:n_ch_down-1),r_ch_down*cos(phi_ch_down),x_tooth(n_ch_down:end)];
   y_tooth= [y_tooth(1:n_ch_down-1),r_ch_down*sin(phi_ch_down),y_tooth(n_ch_down:end)];
end
function V= Mult_homoMat(M,V)
% Matrix multiplication on [4x4] with Vextor area [nx3]
   V= (M * [V,ones(size(V(:,1)))]')';
   V= V(:,1:3);
end
function [p_out,intersec]= Cylinder2line(r_e, p_in, vec)
   % Calculate the insection of a line given by a point of the line
   % and the direction of the line with
   % a cylinder with r_e around the Y-axis
   a_x= p_in(1); a_z= p_in(3);
   n_x= vec(1);  n_z= vec(3);
   tmp= (- a_x^2*n_z^2 + 2*a_x*a_z*n_x*n_z - a_z^2*n_x^2 + n_x^2*r_e^2 + n_z^2*r_e^2);
   if tmp < 0
      intersec= false; p_out=[NaN,NaN,NaN];
   else
      intersec= true;
      l= -(a_x*n_x + a_z*n_z + sqrt(tmp))/(n_x^2 + n_z^2);
      p_out= l*vec + p_in;
   end
end
function [x_tooth, y_tooth]= xy(a,s)
% Calculates with a 4 point cubic spline the point p_ch_down
% add the point p_ch_down to the profile
   %[phi_tooth,r_tooth]= cart2pol(x_tooth(n_ch_down-2:n_ch_down+1),y_tooth(n_ch_down-2:n_ch_down+1));
   %phi_ch_down=spline(r_tooth, phi_tooth, r_ch_down);
   x_tooth=a;
   y_tooth=s;
end
