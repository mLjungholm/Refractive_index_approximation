% v1 = [1 0];
% x1 = [0 0];
% 
% x2 = [1 1];
% v2 = [1 1];
% 
% A = [v1' -v2'];
% y = x2' - x1';
% % if det(A) == 0
% %     a = (x2(1) - x1(1))/v1(1);
% %     if x1(2) + a*v1(2) == x2(2)
% %         disp('infinite')
% %     else
% %         disp('no intersection - parallel')
% %     end
% % else
% is = A\y;
% 
% 
% disp(is)
% ip1 = x1 + is(1).*v1;
% disp(ip1)
% ip2 = x2 + is(2).*v2;
% disp(ip2)
% % end
% 
% 
% close all
% c = 1;
% figure(1)
% hold on; grid on
% quiver(x1(1),x1(2),v1(1),v1(2),c,'r')
% quiver(x2(1),x2(2),v2(1),v2(2),c,'b')
% % plot([p(1); p(1)+c*v(1)], [p(2); p(2)+c*v(2)], 'r')
% % plot([x0(1); x0(1)+c*u(1)], [x0r(2); x0(2)+c*u(2)], 'r')
learningVarargin2 (sin (0 : pi/200 : 4*pi), 'color', 'b')


function learningVarargin2 (x, varargin)
plot (x, varargin {:})
end