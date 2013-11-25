SCC13 graphlab task2
Boston Team
GraphLab README

================



Input Datasets

---------------



- LiveJournal: http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz



- Internet topology: http://snap.stanford.edu/data/as-skitter.txt.gz



- Stiffness matrix:

  http://www.cise.ufl.edu/research/sparse/MM/GHS_psdef/bmw7st_1.tar.gz



- YouTube ground-truth communities:

  http://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz



- Amazon product co-purchasing network:

  http://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz



- Friendster:

  http://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz





TASK 0: Ranking graph nodes on LiveJournal social network data

---------------------------------------------------------------



- Goal: Use pagerank to compute network ranking on the LiveJournal data set

    - Part 1: Install GraphLab as instructed in http://graphlab.org/download

    - Use GraphLab's pagerank to compute network ranking



- What to submit?



    - Console output of the pagerank execution



    - output files containing the pagerank vector



    - The graph node with the highest pagerank value



- Useful links:



    - PageRank method: http://en.wikipedia.org/wiki/PageRank



    - Network Ranking:

      http://docs.graphlab.org/graph_analytics.html#graph_analytics_pagerank



- Tips:

    

    - The LiveJournal dataset format is tsv, remember to set the correct file

      format in the command line



    - Remember to save the results by setting the right command line flag





TASK 1: Spectral Clustering 

----------------------------



- Goal: Perform non-trivial analytics by utilizing GraphLab SVD and k-means



    - Part 1a: run GraphLab SVD on an internet network for performing dimensionality

      reduction of the AS map. The result is the top 10 singular value triplets



    - Part 1b: using the singular vectors of the matrix U, deploy GraphLab k-means

      clustering to compute clusters of Internet AS systems. The result is a

      grouping of related Internet AS systems together

  

    - Part 2a: repeat Part 1a, but this time using the LiveJournal dataset

  

    - Part 2b: repeat Part 1b, but this time using the LiveJournal dataset

  

    - Part 3a: repeat Part 1a, but this time using the YouTube dataset

  

    - Part 3b: repeat Part 1b, but this time using the YouTube dataset



    - Part 4a: repeat Part 1a, but this time using the Amazon dataset

  

    - Part 4b: repeat Part 1b, but this time using the Amazon dataset



- What to submit?

    

    - SVD singular value and error estimation for each of the singular values

      Create a file with the "Singular value ..." lines that appear in the 

      output.



    - k-means run output, including the cost per iteration



- Tips:

    - You need to supply rows and columns of the matrix (graph) size. Use

      --rows=1696416 --cols=1696417 for the AS dataset



    - Use --save_vectors=1 to save SVD output vectors



    - Since the network is binary (namely link weights are all one), use

      --binary=true flag



    - You need to fine tune parameters nsv, nv, max_iter, and tol to get good

      accuracy. When the algorithm is finished you get an estimation of the accuracy

      for each singular value triplet. The nsv parameter is the number of required

      singular values (10 in this case). However, nv should be much bigger to allow

      for a good accuracy. You will need to play with different values of nv until you

      get the required accuracy (1e-12)

   

    - For k-means, you first want to parse the SVD output and create one file

      containing the second column for each output file of the U vectors.

      These columns can then be merged into a single input file that will be

      used for k-means.

      If multiple MPI ranks are used, need to first to merge and sort the

      results.





TASK 2: GraphLab Programming

----------------------------



- Goal: Implement the Jacobi algorithm using GraphLab's warp engine. The Jacobi

  method is a simple parallel algorithm for solving a linear system of

  equations.



    - Part 1: Algorithm implementation



    - Part 2: Algorithm verification on 3x3 matrix (mat3x3.txt)



    - Part 3: Run Jacobi on the Stiffness matrix (bmw7st_1.mtx)



    - Part 4: Run Jacobi on the Friendster data set



- Implementation Requirements:



    - The program should use GraphLab's warp engine.



    - A residual calculation should be performed after each iteration. The

      residual is defined as the second norm squared (|| A*x - b||_2)^2,

      where x is the solution at the end of each iteration.



    - You should support the following command line flags:



        - Number of matrix rows: --rows



        - Number of matrix columns: --cols



        - Max. number of iterations: --max_iter



        - Tolerance threshold: --tol



- Tips:



    - It is advised to follow the warp engine tutorial listed in 'Useful links'

      during implementation



    - The easiest format to work with is the matrix market format



    - Note that the header for the original matrix file must be removed to

      allow you to use GraphLab's graph_loader().



    - You must scale the input matrix, such that the main diagonal will be equal

      to the absolute sum of the row entries +1. The vector b should be all

      ones. See Diagonally dominant matrix in 'Useful links'



    - For an example of how to use clopts, see the svd.cpp source code in the

      collaborative_filtering toolkit



    - The Stiffness matrix is symmetric. You will need to create a GraphLab

      graph out of the sparse symmetric matrix



    - You can use IS_POD_TYPE super class to avoid implementing the save() and

      load() operations which are required when transmitting the vertex data

      over the network. For example:



            struct vertex_data : IS_POD_TYPE { 

                int val1; 

                int val2; 

                vertex_data(int val1, int val2) : val1(val1), val2(val2) { }; 

            }



            typedef vertex_data vertex_data_type;



    - Don't forget to edit the CMakeLists.txt file and add your application. See

      toolkits/linear_solvers/CMakeLists.txt



    - The Jacobi method starts with a random starting vector. It is advised to

      initialize it using the drand48() function. Be careful not to set it to zero!



RECOMMENDATIONS:

    - The Application Specialist used GraphLab commitID 734e9b6f9b to compute

      the results.



    - Task 0 and Task 2 can be completed with the GraphLab v2.2 release



    - It is strongly recommended to use GraphLab commitID 734e9b6f9b for Task

      1, as it contains changes that support binary graphs, and displays cost

      metrics to measure progress and performance.

      The source for that commitID is available in graphlab-734e9b6f9b/ but

      you are welcome to download directly from GitHub



      It is possible to use GraphLab v2.2 for Task 1, but not recommended.



    - A template file has been provided for warp_jacobi_template.cpp to help

      you get started, but not that you are NOT REQUIRED to use it.




