function phi = FcnPhiCSTR(xt,xe,ue,AK,K)
    f = cstr(xt+xe,K*xt+ue);
    phi = f - AK*xt;
end
