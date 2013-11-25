#include <vector>
#include <string>
#include <fstream>
#include <graphlab.hpp>
#include <graphlab/warp.hpp>
#include "../collaborative_filtering/eigen_wrapper.hpp"
#include "../collaborative_filtering/types.hpp"
#include "../collaborative_filtering/eigen_serialization.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <graphlab/util/stl_util.hpp>

using namespace graphlab;



#define initX 0  //random seed for later
//int rows = 3, cols = 3;
int max_iter = 10;
double tol = 1e-5;
//double max_err=100;



// Represents a row of the matrix and its diagonal value
struct vertex_data : public graphlab::IS_POD_TYPE{
    double x_i; // current guess
    //int b_i; // coefficient = side of row
    double A_ii; // diagonal value
    double err;
};

// represents an non-diagonal element of a matrix
struct edge_data : public graphlab::IS_POD_TYPE{
    double A_ij; // value of matrix element
    edge_data(double input = 1) : A_ij(input) {};
};



typedef distributed_graph<vertex_data, edge_data> graph_type;

double max(double x, double y) { return x > y ? x : y; }

/*

double jacobi_map(graph_type::edge_type edge, graph_type::vertex_type other, double prev_guess)
{
    return edge.data().A_ij * prev_guess;
}
*/




double edge_val(graph_type::edge_type edge, graph_type::vertex_type other) { return fabs(edge.data().A_ij); }
void add_magnitude(double& a, const double& b) { a += b;}

void jacobi_precondition(graph_type::vertex_type& vertex)
{
    vertex_data& data = vertex.data();
    double row_sum = warp::map_reduce_neighborhood(vertex, OUT_EDGES, edge_val, add_magnitude);
    data.A_ii = row_sum + 1;
}


double jacobi_map(graph_type::edge_type edge, graph_type::vertex_type other, const double prev_guess) {
    return edge.data().A_ij * prev_guess;
}

void jacobi_combine(double& a, const double& b, const double na) { a += b; }

void jacobi_step(graph_type::vertex_type& vertex)
{
    vertex_data& data = vertex.data();
    double res = warp::map_reduce_neighborhood(vertex, OUT_EDGES, data.x_i, jacobi_map, jacobi_combine);
    double x_i_next = (1 - res) / data.A_ii;
    //if (x_i_next != 0)
	    data.err = fabs(x_i_next - data.x_i) / fabs(x_i_next);
   //else
        //data.err = 0;
    data.x_i = x_i_next;
}



/***
 * JACOBI UPDATE FUNCTION
 * Update rule is:
 * x_i = (b_i - \sum_j A_ij * x_j)/A_ii


void jacobi_step(graph_type::vertex_type& vertex)
{
    vertex_data& data = vertex.data();

    //Warp function allows a map-reduce aggregation of the neighborhood of a vertex to be performed
    double res = warp::map_reduce_neighborhood(vertex, OUT_EDGES, data.x_i, jacobi_map);

    //use OUT_EDGES for sum up the rows
    //IN_EDGES:In edges implies that only whose target is the center vertex are processed during gather or scatter
    //OUT_EDGES:Out edges implies that only whose source is the center vertex are processed during gather or scatter
    //ALL_EDGES:All edges implies that all adges adjacent to a the center vertex are processed on gather or scatter. Note that some neighbors may be encountered twice if there is both an in and out edge to that neighbor

    //b_i =1 always
    //double x_i_next = (data.b_i - res) / data.A_ii;
    double x_i_next = (1 - res) / data.A_ii;
    data.err = fabs(x_i_next - data.x_i) / x_i_next;
    data.x_i = x_i_next;
}


void jacobi_step(graph_type::vertex_type& vertex)
{
    vertex_data& data = vertex.data();
    double res = warp::map_reduce_neighborhood(vertex, OUT_EDGES, data.x_i, jacobi_map, jacobi_combine);
    double x_i_next = (data.b_i - res) / data.A_ii;
    data.err = fabs(x_i_next - data.x_i) / x_i_next;
    data.x_i = x_i_next;
}

  */



//void maxErr(graph_type::vertex_type vertex) {
    //max_err = max(max_err, vertex.data().err);
//}



inline bool graph_loader(graph_type& graph, const std::string& filename, const std::string& line) {
    //no need to parse
    // if (boost::algorithm::ends_with(filename))
    //   return true;

    ASSERT_FALSE(line.empty());
    // Determine the role of the data


    // Parse the line
    std::stringstream strm(line);
    graph_type::vertex_id_type source_id(-1), target_id(-1);
    float weight(0);
    strm >> source_id >> target_id;
    source_id--; target_id--;
    //if (source_id >= (uint)rows)
        //logstream(LOG_FATAL)<<"Row number: " << source_id << " sould be < rows " << rows << " [ line: " << line << " ] " << std::endl;
    //if (target_id >= (uint)cols)
        //logstream(LOG_FATAL)<<"Col number: " << target_id << " sould be < cols " << cols << " [ line: " << line << " ] " << std::endl;
    strm >> weight;

    //the diagonal, which is the vertex value
    if (source_id == target_id){
        vertex_data data;
        data.A_ii = weight; //the diagonal value
        data.x_i = 1.0;
        graph.add_vertex(source_id, data);
    }

    // Create an edge and add it to the graph, for the non diagonals
    else graph.add_edge(source_id, target_id, edge_data(weight));
    return true; // successful load
} // end of graph_loader




struct pagerank_writer{
    std::string save_vertex(graph_type::vertex_type v) {
        char c[128];
        //sprintf(c, "%u\t%f\n", v.id(), v.data().x_i);
        sprintf(c, "%f\n", v.data().x_i);
        return c;
    }
    std::string save_edge(const graph_type::edge_type& edge) {return "";}
};

double err_value(const graph_type::vertex_type& vertex) { return vertex.data().err * vertex.data().err; }



int main(int argc, char** argv)
{
 	graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

    //---------------------------------- input options
 	const std::string description =
    "Solve a linear system using Jacobi method";
    graphlab::command_line_options clopts(description);
	std::string input_dir, output_dir;
    clopts.attach_option("matrix", input_dir,"The directory containing the matrix file");
    clopts.add_positional("matrix");
    //clopts.attach_option("initial_vector", vecfile,"optional initial vector");
    //clopts.attach_option("debug", debug, "Display debug output.");
    clopts.attach_option("max_iter", max_iter, "max iterations");
    clopts.attach_option("tol", tol, "convergence threshold");
    //clopts.attach_option("rows", rows, "number of rows");
    //clopts.attach_option("cols", cols, "number of cols");
    if(!clopts.parse(argc, argv) || input_dir == "") {
        std::cout << "Error in parsing command line arguments." << std::endl;
        clopts.print_description();
        return EXIT_FAILURE;
    }
    //if (rows <= 0 || cols <= 0 || rows != cols)
        //logstream(LOG_FATAL)<<"Please specify number of rows/cols of the input matrix" << std::endl;



    //---------------------------------- calling graph_loader
	dc.cout() << "Loading graph." << std::endl;
    graphlab::timer timer;

    graph_type graph(dc, clopts);				//graph_type
    graph.load(input_dir, graph_loader);     //loading data

    //pgraph = &graph;
    dc.cout() << "Loading graph. Finished in " << timer.current_time() << std::endl;
    dc.cout() << "Finalizing graph." << std::endl;
    timer.start();
    graph.finalize();  //done loading graph
    dc.cout() << "Finalizing graph. Finished in " << timer.current_time() << std::endl;



    //------------------------------------ print graph info
    dc.cout()
    << "========== Graph statistics on proc " << dc.procid() << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replicas() //should be non?
    << "\n Replica to vertex ratio: "
    << float(graph.num_replicas())/graph.num_vertices() //should be 0?
    << "\n --------------------------------------------"
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: "
    << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    //<< "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: "
    << float(graph.num_local_edges())/graph.num_edges()
    << std::endl;


    //----------------------------------
    //

    //warp_engine engine(dc, graph);


    bool converged = false;
    double max_err;
    int k = 0;
    dc.cout() << "Starting Jacobi precondition\n";
   //warp::parfor_all_vertices(graph, jacobi_precondition);
    do {
        dc.cout() << "Starting jacobi iteration " << k << "\n";
        //This Warp function provides a simple parallel for loop over all vertices in the graph, or in a given set of vertices
        warp::parfor_all_vertices(graph, jacobi_step);
        max_err = graph.map_reduce_vertices<double>(err_value);
        max_err = sqrt(max_err);
        //warp::parfor_all_vertices(graph, maxErr);
        dc.cout() << "Error for iteration: " << max_err << std::endl;

        if (fabs(max_err) < 0.000000001) {
            converged = true;
        }
        k++;
    } while (!converged);

    dc.cout() << "Finished computing in " << k << " iterations error: " << max_err << "\n";


     //---------------------------------- finished
     double runtime = timer.current_time();
     dc.cout() << "Jacobi finished in " << runtime << std::endl;
     //dc.cout() << "\t Updates: " << engine.num_updates()<< std::endl;
     //dc.cout() << "Solution converged to residual: " << ret.toDouble() << std::endl;


     dc.cout() << "----------------------------------------------------------"
     << std::endl
     << "Final Runtime (seconds):   " << runtime
     << std::endl
     //<< "Updates executed: " << engine.num_updates()	<< std::endl
     //<< "Update Rate (updates/second): "
     //<< engine.num_updates()	/ (double)runtime << std::endl;
     ;


     //----------------------------------save to file


    //graph.save("gl", pagerank_writer().save_vertex, false, true, false, 1);
    graph.save("gl", pagerank_writer(), false, true, false, 1);


    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;

}
