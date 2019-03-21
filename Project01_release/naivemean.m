%     for i = 1:r.N
%         cluster_i = data(cluster_vector==i,:);
%         r.Centers(1,i) = mean(cluster_i(:,1));
%         r.Centers(2,i) = mean(cluster_i(:,2));
%         r.Color(i) = max(cluster_i(:,3))~=0;
%         [m,~] = size(cluster_i);
%         if m == 1
%             r.Diameter(i) = 0;
%         else     
%             r.Diameter(i) = max(pdist(cluster_i(:,1:2)));
%         end
%     end