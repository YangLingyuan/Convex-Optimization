function projected = func_projection(z,A,b)
Index_active=[0,0];
num_active=0;
%by default, no constraint is violated
projected=z;
for k=1:4
    %if one constraint is violated
    if(A(k,:)*z+b(k)<0)
        num_active=num_active+1;
        Index_active(num_active)=k;
        
        %find the food of perpendicular
        Xp=(A(k,2)*A(k,2)*z(1)-A(k,2)*A(k,1)*z(2)-A(k,1)*b(k))/(A(k,2)*A(k,2)+A(k,1)*A(k,1));
        Yp=(-A(k,2)*A(k,1)*z(1)+A(k,1)*A(k,1)*z(2)-A(k,2)*b(k))/(A(k,2)*A(k,2)+A(k,1)*A(k,1));
        Zp=[Xp,Yp]';
        
        %check if the foot is feasible
        if((A(1,:)*Zp+b(1)>=0)&&(A(2,:)*Zp+b(2)>=0)&&(A(3,:)*Zp+b(3)>=0)&&(A(4,:)*Zp+b(4)>=0))
            %if feet of the perpendicular is within the feasible set,
            %it should be the optimal solution choose it as the projected point
            projected=Zp;
            break;
        end
    
    
    %if z violate 2 constraints and both foot of the
    %perpendicular is not within the feasible set, choose the
    %intersection point as optimal
    if(num_active==2)
        
        z1=(-sign(A(Index_active(2),2))*b(Index_active(2))--sign(A(Index_active(1),2))*b(Index_active(1)))/(A(Index_active(2),1)/A(Index_active(2),2)-A(Index_active(1),1)/A(Index_active(1),2));
        z2=-A(Index_active(1),1)/A(Index_active(1),2)*z1-sign(A(Index_active(1),2))*b(Index_active(1));
        projected=[z1,z2]';
        break;
    end
    
    %in the end if only one constraint is violated and its foot is not
    %feasible
    if(k==4 && num_active==1)
        for j=1:4
            if(A(j,:)*Zp+b(j)<0)
                z1=(-sign(A(j,2))*b(j)--sign(A(Index_active(1),2))*b(Index_active(1)))/(A(j,1)/A(j,2)-A(Index_active(1),1)/A(Index_active(1),2));
                z2=-A(Index_active(1),1)/A(Index_active(1),2)*z1+-sign(A(Index_active(1),2))*b(Index_active(1));
                projected=[z1,z2]';
                break;
            end
        end
    end
    end
end

end