function [iip,Kip, iim,Kim, jjp,Kjp, jjm,Kjm,I,Ip,Im,J,Jp,Jm,Mis,Pls]=getiijj(edge,tri_for_e,edges_for_tri,orientation_nor_tri,orientation_tang_tri,hex_for_e,hex_A_section,edges_for_hex,tri_for_hex)


if(orientation_nor_tri(1,edge)>0)
    I=tri_for_e(2,edge);
    J=tri_for_e(1,edge);
else
    I=tri_for_e(1,edge);
    J=tri_for_e(2,edge);
end


    if(orientation_tang_tri(1,edge)<0)
        Mis=hex_for_e(2,edge);
        Pls=hex_for_e(1,edge);
    else
        Mis=hex_for_e(1,edge);
        Pls=hex_for_e(2,edge);
    end


iip=setdiff(edges_for_tri(:,I),edges_for_hex(:,Mis));
iim=setdiff(edges_for_tri(:,I),edges_for_hex(:,Pls));
jjp=setdiff(edges_for_tri(:,J),edges_for_hex(:,Mis));
jjm=setdiff(edges_for_tri(:,J),edges_for_hex(:,Pls));

Ip=tri_for_e(:,iip); Ip(Ip==I)=[];
Im=tri_for_e(:,iim); Im(Im==I)=[];
Jp=tri_for_e(:,jjp); Jp(Jp==J)=[];
Jm=tri_for_e(:,jjm); Jm(Jm==J)=[];

Kip=hex_A_section(find(tri_for_hex(:,Pls)==I),Pls);
Kim=hex_A_section(find(tri_for_hex(:,Mis)==I),Mis);
Kjp=hex_A_section(find(tri_for_hex(:,Pls)==J),Pls);
Kjm=hex_A_section(find(tri_for_hex(:,Mis)==J),Mis);

end
