function node_index3D=chuyennutvaophantu3D(spt3D,luoi_phantu3D,luoi3D)
node_index3D=zeros(8,spt3D);
for i=1:spt3D
    
    [R,S,T]=ind2sub(size(luoi_phantu3D),find(luoi_phantu3D==i));
    
    node_index3D(1:8,i)=[luoi3D((R+1),S,T),luoi3D((R+1),S+1,T),luoi3D(R,S+1,T),luoi3D(R,S,T),...
            luoi3D((R+1),S,T+1),luoi3D((R+1),S+1,T+1),luoi3D(R,S+1,T+1),luoi3D(R,S,T+1)];
    
end