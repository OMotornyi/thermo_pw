import settings;
defaultpen(1.6);
import three;
currentprojection=orthographic(   5.000,   2.000,   1.000);
size(15cm);
size3(40cm,40cm,40cm);

triple    F1=(  1.000000,  0.000000,  0.500000);
triple    F2=(  1.000000,  0.500000, -0.000000);
triple    F3=(  1.000000,  0.000000, -0.500000);
triple    F4=(  1.000000, -0.500000, -0.000000);
triple    F5=(  0.000000,  1.000000,  0.500000);
triple    F6=( -0.500000,  1.000000, -0.000000);
triple    F7=(  0.000000,  1.000000, -0.500000);
triple    F8=(  0.500000,  1.000000, -0.000000);
triple    F9=( -1.000000,  0.000000,  0.500000);
triple   F10=( -1.000000, -0.500000, -0.000000);
triple   F11=( -1.000000,  0.000000, -0.500000);
triple   F12=( -1.000000,  0.500000, -0.000000);
triple   F13=( -0.000000, -1.000000,  0.500000);
triple   F14=(  0.500000, -1.000000, -0.000000);
triple   F15=( -0.000000, -1.000000, -0.500000);
triple   F16=( -0.500000, -1.000000, -0.000000);
triple   F17=(  0.500000, -0.000000,  1.000000);
triple   F18=(  0.000000,  0.500000,  1.000000);
triple   F19=( -0.500000, -0.000000,  1.000000);
triple   F20=( -0.000000, -0.500000,  1.000000);
triple   F21=(  0.500000, -0.000000, -1.000000);
triple   F22=(  0.000000,  0.500000, -1.000000);
triple   F23=( -0.500000, -0.000000, -1.000000);
triple   F24=( -0.000000, -0.500000, -1.000000);

path3 yp=F1--F2--F3--F4--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F5--F6--F7--F8--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F17--F18--F19--F20--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F1--F2--F8--F5--F18--F17--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F1--F17--F20--F13--F14--F4--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F3--F2--F8--F7--F22--F21--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);

path3 yp=F3--F4--F14--F15--F24--F21--cycle;
surface y=surface(yp);
draw(y,blue+opacity(0.7));
draw(yp);
triple   PP1=(  0.000000,  0.000000,  0.000000);
label(scale(1.6)*"$\Gamma$",PP1,S,red);
triple   PP2=(  0.000000,  1.000000,  0.000000);
label(scale(1.6)*"X",PP2,NE,red);
triple   PP3=(  0.500000,  1.000000, -0.000000);
label(scale(1.6)*"W",PP3,SE,red);
triple   PP4=(  1.000000,  1.000000,  0.000000);
label(scale(1.6)*"X",PP4,SE,red);
triple   PP5=(  0.750000,  0.750000, -0.000000);
label(scale(1.6)*"K",PP5,S,red);
triple   PP6=(  0.000000,  0.000000,  0.000000);
label(scale(1.6)*"$\Gamma$",PP6,S,red);
triple   PP7=(  0.500000,  0.500000,  0.500000);
label(scale(1.6)*"L",PP7,N,red);
triple      MM1=(  0.000000,  1.000000,  0.000000);
triple      MM2=(  0.000000,  0.000000,  0.000000);
draw(MM1--MM2,red+dashed);
triple      MM3=(  0.500000,  1.000000, -0.000000);
triple      MM4=(  0.000000,  1.000000,  0.000000);
draw(MM3--MM4,red);
triple      MM5=(  1.000000,  1.000000,  0.000000);
triple      MM6=(  0.500000,  1.000000, -0.000000);
draw(MM5--MM6,red);
triple      MM7=(  0.750000,  0.750000, -0.000000);
triple      MM8=(  1.000000,  1.000000,  0.000000);
draw(MM7--MM8,red);
triple      MM9=(  0.000000,  0.000000,  0.000000);
triple     MM10=(  0.750000,  0.750000, -0.000000);
draw(MM9--MM10,red+dashed);
triple     MM11=(  0.500000,  0.500000,  0.500000);
triple     MM12=(  0.000000,  0.000000,  0.000000);
draw(MM11--MM12,red+dashed);
triple G=(0.0,0.0,0.0);
triple M1=(  1.000000,0.0,0.0);
triple M2=(0.0,  1.000000,0.0);
triple M3=(0.0,0.0,  1.000000);
draw(G--M1,dotted);
draw(M1--M1+(  0.550000,0.0,0.0),Arrow3);
draw(G--M2,dotted);
draw(M2--M2+(0.0,  0.450000,0.0),Arrow3);
draw(G--M3,dotted);
draw(M3--M3+(0.0,0.0,  0.550000),Arrow3);

label(scale(1.9)*"$k_x$",M1+(  0.550000,0.0,0.0),NW);
label(scale(1.9)*"$k_y$",M2+(0.0,  0.450000,0.0),N);
label(scale(1.9)*"$k_z$",M3+(0.0,0.02,  0.550000),SE);

