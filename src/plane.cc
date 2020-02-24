#include <unordered_map>
#include <queue>
#include <set>

#include <igl/readOBJ.h>

#define GRID_RESOLUTION 1000
#define ANGLE_THRES 5.0

typedef std::pair<int, std::pair<int, int> > VertexHash;
void ComputeSurfaceNormals(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd* pFN) {
	*pFN = Eigen::MatrixXd(F.rows(), 3);
	Eigen::MatrixXd& FN = *pFN;
	for (int i = 0; i < F.rows(); ++i) {
		Eigen::Vector3d v0 = V.row(F(i, 0));
		Eigen::Vector3d v1 = V.row(F(i, 1));
		Eigen::Vector3d v2 = V.row(F(i, 2));
		Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
		double l = n.norm();
		if (l > 0) {
			n = n / l;
		}
		FN.row(i) = n;
	}
}

void NormalizeVertices(Eigen::MatrixXd& V) {
	double min_corner[] = {1e30, 1e30, 1e30};
	double max_corner[] = {-1e30, -1e30, -1e30};
	for (int i = 0; i < V.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			double v = V(i, j);
			min_corner[j] = std::min(min_corner[j], v);
			max_corner[j] = std::max(max_corner[j], v);
		}
	}
	for (int i = 0; i < V.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			double v = V(i, j);
			v = (v - min_corner[j]) / (max_corner[j] - min_corner[j] + 1e-6);
			V(i, j) = v;
		}
	}
}

void VisualizeSurfaceColor(const char* filename,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& FN) {
	std::ofstream os(filename);
	for (int i = 0; i < F.rows(); ++i) {
		Eigen::Vector3d v0 = V.row(F(i, 0));
		Eigen::Vector3d v1 = V.row(F(i, 1));
		Eigen::Vector3d v2 = V.row(F(i, 2));
		Eigen::Vector3d n = FN.row(i);		
		n = n * 0.49 + Eigen::Vector3d(0.5, 0.5, 0.5);
		os << "v " << v0[0] << " " << v0[1] << " " << v0[2] << " "
		   << n[0] << " " << n[1] << " " << n[2] << "\n";
		os << "v " << v1[0] << " " << v1[1] << " " << v1[2] << " "
		   << n[0] << " " << n[1] << " " << n[2] << "\n";
		os << "v " << v2[0] << " " << v2[1] << " " << v2[2] << " "
		   << n[0] << " " << n[1] << " " << n[2] << "\n";
	}

	for (int i = 0; i < F.rows(); ++i) {
		os << "f " << i * 3 + 1 << " " << i * 3 + 2 << " " << i * 3 + 3 << "\n";
	}

	os.close();
}

VertexHash EncodeVertex(const Eigen::Vector3d& v) {
	int x = v[0] * GRID_RESOLUTION;
	int y = v[1] * GRID_RESOLUTION;
	int z = v[2] * GRID_RESOLUTION;
	return std::make_pair(x, std::make_pair(y, z));
}

long long Key(const VertexHash& h) {
	return ((long long)h.first * GRID_RESOLUTION + (long long)h.second.first)
		* GRID_RESOLUTION + (long long)h.second.second;
}

void CreateNeighborHash(const VertexHash& h, std::vector<VertexHash>* h_array) {
	for (int i = -1; i < 2; ++i) {
		for (int j = -1; j < 2; ++j) {
			for (int k = -1; k < 2; ++k) {
				VertexHash nh = h;
				nh.first += i;
				nh.second.first += j;
				nh.second.second += k;
				h_array->push_back(nh);
			}
		}
	}
}

int main(int argc, char** argv)
{
	if (argc < 2) {
		printf("./plane input.obj output.txt [visual.obj]\n");
		return 0;
	}
	const char* input_file = argv[1];
	const char* output_file = argv[2];
	const char* visual_file = (argc > 3) ? argv[3] : 0;

	Eigen::MatrixXd V_orig;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd FN;
	igl::readOBJ(input_file, V, F);

	if (argc > 3)
		V_orig = V;
	NormalizeVertices(V);
	ComputeSurfaceNormals(V, F, &FN);

	//VisualizeSurfaceNormal(output_file, V, F, FN);

	printf("hash faces with vertex location!\n");
	std::unordered_map<long long, std::vector<int> > vertex_to_faces;

	for (int i = 0; i < F.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			Eigen::Vector3d v = V.row(F(i, j));
			VertexHash h = EncodeVertex(v);
			long long k = Key(h);
			auto it = vertex_to_faces.find(k);
			if (it == vertex_to_faces.end()) {
				std::vector<int> faces;
				faces.push_back(i);
				vertex_to_faces[k] = faces;
			} else {
				it->second.push_back(i);
			}
		}
	}

	std::vector<int> visited(F.rows(), 0);
	int num_group = 1;
	for (int i = 0; i < F.rows(); ++i) {
		if (visited[i] != 0)
			continue;
		Eigen::Vector3d n = FN.row(i);
		if (n.norm() == 0)
			continue;
		visited[i] = num_group;
		std::queue<int> faces_queue;
		faces_queue.push(i);
		while (!faces_queue.empty()) {
			int face_id = faces_queue.front();
			Eigen::Vector3d n = FN.row(face_id);
			faces_queue.pop();

			for (int j = 0; j < 3; ++j) {
				Eigen::Vector3d v = V.row(F(face_id, j));
				VertexHash h = EncodeVertex(v);
				std::vector<VertexHash> h_array;
				CreateNeighborHash(h, &h_array);

				for (auto& nh : h_array) {
					long long k = Key(nh);
					auto it = vertex_to_faces.find(k);
					if (it == vertex_to_faces.end())
						continue;
					for (auto& nf : it->second) {
						if (visited[nf])
							continue;
						Eigen::Vector3d nn = FN.row(nf);
						if (std::acos(std::abs(nn.dot(n))) / 3.141592654 * 180.0 < ANGLE_THRES) {
							visited[nf] = num_group;
							faces_queue.push(nf);
						}
					}
				}
			}
		}
		num_group += 1;
	}

	std::ofstream os(output_file);
	for (int i = 0; i < visited.size() - 1; ++i) {
		os << visited[i] << " ";
	}
	os << visited.back() << "\n";
	os.close();

	if (visual_file) {

		std::vector<Eigen::Vector3d> color_groups(num_group);
		for (int i = 0; i < num_group; ++i) {
			color_groups[i] = Eigen::Vector3d(rand() / (double)RAND_MAX * 2.0 - 1.0,
				rand() / (double)RAND_MAX * 2.0 - 1.0,
				rand() / (double)RAND_MAX * 2.0 - 1.0);
		}
		color_groups[0] = Eigen::Vector3d(0, 0, 0);

		Eigen::MatrixXd C(visited.size(), 3);
		for (int i = 0; i < visited.size(); ++i) {
			C.row(i) = color_groups[visited[i]];
		}

		VisualizeSurfaceColor(visual_file, V_orig, F, C);
	}
	return 0;
}

