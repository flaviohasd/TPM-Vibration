function a = aceleracaotpm(estados,K,C,F,me,d1,d2,K1,K2)
    a(1) = estados(2);
    a(2) = (1/me)*(-K*estados(1) -C*estados(2) + F + K1*d1 + K2*d2);
end