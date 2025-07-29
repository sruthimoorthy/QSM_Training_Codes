function plot_pointcloud_with_qsm(P, cylinder, fig, nf, alp, Bal, Sub, Ind)
% plot_pointcloud_with_qsm  Overlay TLS point‐cloud P with QSM cylinders.
%
%   P         : N×3 (or N×2) point‐cloud
%   cylinder  : struct array or numeric array of cylinder fits
%   fig       : figure number

%   nf        : max facets in the thickest cylinder (min 4)
%   alp       : patch transparency (0–1)
%   Bal, Sub  : cover‐set selection (as in point_cloud_plotting)
%   Ind       : indices of cylinders to plot (optional)

  %— handle defaults for point size & transparency —%
  
  if nargin<4
    error('Need nf (number of facets).')
  end
  if nargin<5 || isempty(alp), alp = 0.6; end

  %— select point indices I based on Bal/Sub —%
  if nargin < 6
    I = 1:size(P,1);
  elseif nargin==7
    I = vertcat(Bal{:});
  else
    if iscell(Sub)
      tmp = vertcat(Sub{:});
      tmp = vertcat(tmp{:});
      I = vertcat(Bal{tmp});
    else
      I = vertcat(Bal{Sub});
    end
  end

  %— set up figure and plot the point cloud —%
figure(fig); clf; hold on
  pcAlpha = 0.2;               % adjust transparency from 0 (invisible) to 1 (opaque)
  pcColor = [0.5 0.5 0.5];     % light gray
  ms = 0.8;                     % marker size for points 
  if size(P,2)==3
    % use scatter3 for alpha support
    scatter3(P(I,1), P(I,2), P(I,3), ms, ...
             'MarkerFaceColor', pcColor, ...
             'MarkerEdgeColor', pcColor, ...
             'MarkerFaceAlpha', pcAlpha, ...
             'MarkerEdgeAlpha', pcAlpha);
  else
    scatter(P(I,1), P(I,2), ms, ...
            'MarkerFaceColor', pcColor, ...
            'MarkerEdgeColor', pcColor, ...
            'MarkerFaceAlpha', pcAlpha, ...
            'MarkerEdgeAlpha', pcAlpha);
  end
  axis equal

  %— BEGIN your original plot_cylinder_model code —%

  if isstruct(cylinder)
      Rad = cylinder.radius;
      Len = cylinder.length;
      Sta = cylinder.start;
      Axe = cylinder.axis;
      BOrd = cylinder.BranchOrder;
  else
      Rad = cylinder(:,1);
      Len = cylinder(:,2);
      Sta = cylinder(:,3:5);
      Sta = mat_vec_subtraction(Sta,Sta(1,:));
      Axe = cylinder(:,6:8);
      BOrd = cylinder(:,9);
  end

  if nargin == 9
      Rad = Rad(Ind);
      Len = Len(Ind);
      Sta = Sta(Ind,:);
      Axe = Axe(Ind,:);
      BOrd = BOrd(Ind);
  end

  nc = size(Rad,1);

  col = [
      0.00  0.00  1.00
      0.00  0.50  0.00
      1.00  0.00  0.00
      0.00  0.75  0.75
      0.75  0.00  0.75
      0.75  0.75  0.00
      0.25  0.25  0.25
      0.75  0.25  0.25
      0.95  0.95  0.00
      0.25  0.25  0.75
      0.75  0.75  0.75
      0.00  1.00  0.00
      0.76  0.57  0.17
      0.54  0.63  0.22
      0.34  0.57  0.92
      1.00  0.10  0.60
      0.88  0.75  0.73
      0.10  0.49  0.47
      0.66  0.34  0.65
      0.99  0.41  0.23];

  N = double(max(BOrd)+1);
  if N <= 20
      col = col(1:N,:);
  else
      m = ceil(N/20);
      col = repmat(col,[m,1]);
      col = col(1:N,:);
  end

  Cir = cell(nf,2);
  for i = 4:nf
      B = [cos((1/i:1/i:1)*2*pi)' sin((1/i:1/i:1)*2*pi)' zeros(i,1)];
      T = [cos((1/i:1/i:1)*2*pi)' sin((1/i:1/i:1)*2*pi)' ones(i,1)];
      Cir{i,1} = [B; T];
      Cir{i,2} = [(1:1:i)' (i+1:1:2*i)' [(i+2:1:2*i)'; i+1] [(2:1:i)'; 1]];
  end

  Vert = zeros(2*nc*(nf+1),3);
  Facets = zeros(nc*(nf+1),4);
  fvd = zeros(nc*(nf+1),3);
  t = 1;
  f = 1;

  for i = 1:nc
      n = ceil(sqrt(Rad(i)/Rad(1))*nf);
      n = min(n,nf);
      n = max(n,4);
      C = Cir{n,1};
      % Scale
      C(:,1:2) = Rad(i)*C(:,1:2);
      C(n+1:end,3) = Len(i)*C(n+1:end,3);
      % Rotate
      ang = real(acos(Axe(i,3)));
      Axis = cross([0 0 1]',Axe(i,:)');
      Rot = rotation_matrix(Axis,ang);
      C = (Rot*C')';
      % Translate
      C = mat_vec_subtraction(C,-Sta(i,:));
      Vert(t:t+2*n-1,:) = C;
      Facets(f:f+n-1,:) = Cir{n,2}+t-1;
      fvd(f:f+n-1,:) = repmat(col(BOrd(i)+1,:),[n 1]);
      t = t+2*n;
      f = f+n;
  end
  t = t-1;
  f = f-1;
  Vert = Vert(1:t,:);
  Facets = Facets(1:f,:);
  fvd = fvd(1:f,:);

  % plot the combined cylinder mesh
  patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
  alpha(alp)
  axis equal
  grid on
  view(-37.5,30)

  %— END original code —%

end

function R = rotation_matrix(u, theta)
  u = u ./ norm(u);
  c = cos(theta); s = sin(theta); C = 1-c;
  ux = u(1); uy = u(2); uz = u(3);
  R = [ ux*ux*C + c,    ux*uy*C - uz*s, ux*uz*C + uy*s;
        uy*ux*C + uz*s, uy*uy*C + c,    uy*uz*C - ux*s;
        uz*ux*C - uy*s, uz*uy*C + ux*s, uz*uz*C + c ];
end
