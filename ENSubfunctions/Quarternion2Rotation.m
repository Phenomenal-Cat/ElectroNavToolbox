% Calculate rotation matrix from quaternions in Nifti header

function R = Quarternion2Rotation(b,c,d)

a = sqrt(1-b^2-c^2-d^2);

R = [   a^2+b^2-c^2-d^2, 	2*(b*c - a*d),      2*(b*d + a*c);...
        2*(b*c + a*d),      a^2+c^2-b^2-d^2,    2*(c*d - a*b);...
        2*(b*d - a*c),      2*(c*d + a*b),      a^2+d^2-b^2-c^2];

end