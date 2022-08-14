function a = aceleracao(estados,K,C,F,me)
    a(1) = estados(2);
    a(2) = -K*estados(1)/me -C*estados(2)/me + F/me;
end
