using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.Drawing;
using System.Diagnostics;
using WindowsFormsApplication1;


namespace TSP
{

    class ProblemAndSolver
    {  
        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;

		private CityMap cityMap;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            // Cleverly come up with an initial BSSF
            // Start at City[0]
            Stopwatch timer = new Stopwatch();
            timer.Start();
            ArrayList initial = new ArrayList();
            initial.Add(Cities[0]);

            List<City> temp = new List<City>();
            for (int i = 0; i < Cities.Length; i++)                             // Time O(n)
                temp.Add(Cities[i]);                                            // Space O(n)
            
            HeapQue initialQue = new HeapQue(temp);                             // Space O(n) Time O(n)
            initialQue.decreaseKey(0, 0);                                       // Time O(1)
            for (int i = 1; i < Cities.Length; i++)
                initialQue.decreaseKey(i, Cities[0].costToGetTo(Cities[i]));    // Time O(logn)
            
            while(!initialQue.isEmpty())
                initial.Add(Cities[initialQue.deleteMin()]);

            bssf = new TSPSolution(initial);

            // Begin the Branch and Bound
			
            // Produce the Initial Reduced Cost Matrix
            double[,] rcm = new double[Cities.Length, Cities.Length];           // Space O(n^2)
            double lowerBound = 0;
            for(int y = 0; y < Cities.Length; y++)                              // Time O(n^2)
            {
                double shortestDistance = Double.PositiveInfinity;
                // initialize the matrix with the distance to each other city x from city y, keep track of the closet city
                for(int x = 0; x < Cities.Length; x++)                          // Time O(n)
                {
                    if (x == y)
                        rcm[x, y] = Double.PositiveInfinity;
                    else
                        rcm[x, y] = Cities[x].costToGetTo(Cities[y]);
                    if (rcm[x, y] < shortestDistance)
                        shortestDistance = rcm[x, y];
                }

                // add the distance to the closest city to the lower bound and subtract the distance to the closest city from the distance to every other city
                lowerBound += shortestDistance;
                for(int x = 0; x < Cities.Length; x++)                          // Time O(n)
                {
                    rcm[x, y] = rcm[x, y] - shortestDistance;
                }
            }
            // make sure there is at least one zero in each column, or that each column is all infinity
            for(int x = 0; x < Cities.Length; x++)                              // Time O(n^2)
            {
                double shortestDistance = Double.PositiveInfinity;
                for(int y = 0; y < Cities.Length; y++)                          // Time O(n)
                {
                    if (rcm[x, y] < shortestDistance)
                        shortestDistance = rcm[x, y];
                }
                if(shortestDistance != 0 && shortestDistance != Double.PositiveInfinity)
                {
                    // if the shortest distance in the column is not zero and it is not infinity, add the smallest value to the lower bound and decrement every other value in the column
                    lowerBound += shortestDistance;
                    for(int y = 0; y < Cities.Length; y++)                      // Time O(n)
                    {
                        rcm[x, y] = rcm[x, y] - shortestDistance;
                    }
                }
            }

			// Priority Que
			PriorityQueue<BBState> que = new PriorityQueue<BBState>();

			// Find the Arnell Constant
			double arnellConstant = 0;
			int numPaths = 0;
			for(int y = 0; y < Cities.Length; y++)
			{
				for(int x = 0; x < Cities.Length; x++)
				{
					double cost = Cities[x].costToGetTo(Cities[y]);
					if (cost < Double.PositiveInfinity)
					{
						arnellConstant += cost;
						numPaths++;
					}
				}
			}
			arnellConstant = arnellConstant / numPaths;
			BBState.arnellConstant = arnellConstant;
			Console.WriteLine("This is not broken");



			// all counting variables
			int count = 0;
			int stateCount = 0;
			int trimCount = 0;
			int storedCount = 0;
			String bestTime = timer.Elapsed.ToString();

			// add states to the Priority Que representing starting at each city
			for(int i = 0; i < Cities.Length; i++)                              // Time O(n)
			{
				List<int> newRoute = new List<int>();                           // Space O(n)
				newRoute.Add(i);
				double[,] newRCM = (double[,])rcm.Clone();						// Space O(n^2)
				que.Enqueue(new BBState(newRCM, lowerBound, newRoute));			// Space O(n)
				stateCount++;
			}
			
			// while there is still time left, keep looping
			while(timer.Elapsed.TotalMilliseconds < time_limit && que.Count() > 0)  // Space O(n!n^3) Time O(n!n^3)  
			{
				// Count the max number of states stored in the queue at one time
				if (storedCount < que.Count())
					storedCount = que.Count();

				// Get the next state, if its lower bound is greater than the BSSF prune it
				BBState state = que.Dequeue();
				if (state.lowerBound > bssf.costOfRoute())
				{
					trimCount++;
					continue;
				}
				// Find all cities that are yet to be visited
				HashSet<int> unvisitedCities = new HashSet<int>();
				for(int i = 0; i < Cities.Length; i++)							// Time O(n)
				{
					unvisitedCities.Add(i);										// Space O(n)
				}
				for(int i = 0; i < state.route.Count; i++)						// Time O(n)
				{
					unvisitedCities.Remove((int)state.route[i]);
				}
				foreach(int nextCity in unvisitedCities)						// Time O(n^3) Space O(n^3)
				{
					int startingCity = (int)state.route[state.route.Count - 1];
					BBState nextState = GenerateNextState(rcm, state.route, nextCity, state.lowerBound);  // Space O(n^2) Time O(n^2)
					stateCount++;
					if (nextState.lowerBound > bssf.costOfRoute())
					{
						trimCount++;
						continue;
					}
					// check if the route is complete, if it is and better than BSSF, set as BSSF
					if (nextState.route.Count == Cities.Length)
					{
						ArrayList cityRoute = new ArrayList();
						for(int i = 0; i < nextState.route.Count; i++)			// Time O(n)
						{
							cityRoute.Add(Cities[nextState.route[i]]);			// Space O(n)
						}
						TSPSolution newSolution = new TSPSolution(cityRoute);
						if (newSolution.costOfRoute() < bssf.costOfRoute())
						{
							bssf = newSolution;
							bestTime = timer.Elapsed.ToString();
							count++;
						}
					}
					else // if the route is not complete, add it to the priority queue
					{
						que.Enqueue(nextState);
					}
				}
			}
			timer.Stop();

			// Print count variables
			Console.WriteLine("Run Time: " + bestTime + " Max stored#: " + storedCount + " Total Created: " + stateCount + " Pruned#: " + (trimCount + que.Count()) + " Arnell's Const: " + BBState.arnellConstant);

			// pass array of results to GUI
            string[] results = new string[3];
            results[COST] = costOfBssf().ToString();    
            results[TIME] = bestTime;
            results[COUNT] = count + "";

            return results;
        }

		// create the next state representing going from the starting city to the next city
		private BBState GenerateNextState(double[,] reducedCostMatrix, List<int> currentRoute, int nextCity, double lowerBound)
		{
			// Time O(n^2) Space O(n^2) 
			// find starting city
			List<int> route = new List<int>();
			for(int i = 0; i < currentRoute.Count; i++)														// Time O(n)
			{
				route.Add(currentRoute[i]);																	// Space O(n)
			}
			int startingCity = (int)route[route.Count - 1];
			// add the next city to the route
			route.Add(nextCity);
			// copy the reduced cost matrix to a new rcm
			double[,] rcm = (double[,])reducedCostMatrix.Clone();											// Space O(n^2)
			// add the value of going from the starting city to next city to the current lower bound
			lowerBound += rcm[startingCity, nextCity];
			// make the inverse value infinity
			rcm[nextCity, startingCity] = Double.PositiveInfinity;

			// make the starting city row and next city column all infinity
			for(int i = 0; i < Cities.Length; i++)															// Time O(n)
			{
				rcm[startingCity, i] = Double.PositiveInfinity;
				rcm[i, nextCity] = Double.PositiveInfinity;
			}

			// rebalance the array, making sure that every row and column has at least one 0 or is all infinity
			for (int x = 0; x < Cities.Length; x++)															// Time O(n^2)
			{
				double shortestDistance = Double.PositiveInfinity;
				for (int y = 0; y < Cities.Length; y++)														// Time O(n)
				{
					if (rcm[x, y] < shortestDistance)
						shortestDistance = rcm[x, y];
				}
				if (shortestDistance != 0 && shortestDistance != Double.PositiveInfinity)
				{
					// if the shortest distance in the column is not zero and it is not infinity, add the smallest value to the lower bound and decrement every other value in the column
					lowerBound += shortestDistance;
					for (int y = 0; y < Cities.Length; y++)													// Time O(n)
					{
						rcm[x, y] = rcm[x, y] - shortestDistance;
					}
				}
			}

			return new BBState(rcm, lowerBound, route);
		}

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
			Console.WriteLine("Running greedy solutions");
			Stopwatch timer = new Stopwatch();
			timer.Start();

			// set of all greedy solutions
			HashSet<TSPSolution> greedySolutions = new HashSet<TSPSolution>();

			// find the greedy solution for the path starting at each city
			for(int i = 0; i < Cities.Length; i++)
			{
				Console.WriteLine("Finding " + i);
				List<int> route = new List<int>();
				// at first all cities are unvisited
				HashSet<int> unvisited = new HashSet<int>();
				for (int j = 0; j < Cities.Length; j++)
				{
					unvisited.Add(j);
				}

				route.Add(i);
				unvisited.Remove(i);

				while(unvisited.Count > 0)
				{
					int closestIndex = -1;
					double closestDistance = Double.PositiveInfinity;

					// find the closest city to the last city in the current route
					foreach(int city in unvisited)
					{
						double distance = Cities[route[route.Count - 1]].costToGetTo(Cities[city]);
						if(distance < closestDistance)
						{
							closestIndex = city;
							closestDistance = distance;
						}
					}

					// add that closest city to the route
					route.Add(closestIndex);
					unvisited.Remove(closestIndex);
				}
				// greedy solution is complete, add it to the set of all solutions
				ArrayList solution = new ArrayList();
				for(int j = 0; j < route.Count; j++)
				{
					solution.Add(Cities[route[j]]);
				}
				TSPSolution tspSolution = new TSPSolution(solution);
				greedySolutions.Add(tspSolution);
				Console.WriteLine("Cost: " + tspSolution.costOfRoute());
			}

			// find the best greedy solution
			TSPSolution greediest = null;
			foreach(TSPSolution solution in greedySolutions)
			{
				if(greediest == null || greediest.costOfRoute() > solution.costOfRoute())
				{
					greediest = solution;
				}
			}
			bssf = greediest;
			timer.Stop();

            string[] results = new string[3];
            results[COST] = greediest.costOfRoute().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = Cities.Length + "";

            return results;
        }

        public string[] fancySolveProblem()
        {
			// Implement a genetic algorithm
			// set global data for Genetic Solution Citizens
			GSCitizen.setCities(Cities);
			GSCitizen.setRandom(rnd);
			int solutionsCount = 0;
			int generationCount = 1;
			// start time
			Stopwatch timer = new Stopwatch();
			timer.Start();
			// Best Citizen So Far
			GSCitizen bcsf = getRandomSolution();
			solutionsCount++;
			String bestTime = timer.Elapsed.ToString();
			// set global data for Genetic Solution
			int populationSize = Cities.Length * 10;
			int groupSize = 5;
			int numNearbyCities = Cities.Length;
			cityMap = new CityMap(Cities, numNearbyCities);
			List<GSCitizen> population = new List<GSCitizen>();

			// find initial greedy solutions
			for(int i = 0; i < populationSize; i++)
			{
				GSCitizen citizen = getGreedySolution(i % Cities.Length);
				if (bcsf.fitness() > citizen.fitness())
				{
					bcsf = citizen;
					bestTime = timer.Elapsed.ToString();
					solutionsCount++;
				}
				population.Add(citizen);
			}

			// Start breeding
			while(timer.Elapsed.TotalMilliseconds < time_limit)
			//while(generationCount < 200)
			{
				// produce the next generation
				List<GSCitizen> nextGeneration = new List<GSCitizen>();
				// shuffle the population
				Shuffle<GSCitizen>(population);
				// take groups of groupSize, breed the fitest two and replace the least fit group members with surviving children
				while(population.Count > 0)
				{
					SpawningPool pool = new SpawningPool();
					for(int i = 0; i < groupSize; i++)
					{
						if (population[0] != null)
						{
							pool.Add(population[0]);
							population.RemoveAt(0);
						}
					}
					GSCitizen[] survivors = pool.spawn();
					for(int i = 0; i < survivors.Length; i++)
					{
						nextGeneration.Add(survivors[i]);
						// if any survivers are more fit than the best citizen so far, make it the new best citizen so far
						if (survivors[i].fitness() < bcsf.fitness())
						{
							bcsf = survivors[i];
							bestTime = timer.Elapsed.ToString();
							solutionsCount++;
						}
					}
				}
				population = nextGeneration;
				generationCount++;
				Console.WriteLine("Finished a generation");
			}

			//bssf = getGreedySolution(0).getSolution();
			bssf = bcsf.getSolution();
			double[] genValues = new double[populationSize];
			for(int i = 0; i < populationSize; i++)
			{
				genValues[i] = population[i].fitness();
			}

            string[] results = new string[3];
            results[COST] = costOfBssf().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = bestTime;
			results[COUNT] = solutionsCount + "";

            return results;
        }

		private GSCitizen getGreedySolution(int startCity)
		{
			List<int> route = new List<int>();
			// at first all cities are unvisited
			HashSet<int> unvisited = new HashSet<int>();
			for (int j = 0; j < Cities.Length; j++)
			{
				unvisited.Add(j);
			}

			route.Add(startCity);
			unvisited.Remove(startCity);

			while (unvisited.Count > 0)
			{
				int next = -1;
				int[] nearby = cityMap.getClosestCities(route[route.Count - 1]);
				bool available = false;
				for (int i = 0; i < nearby.Length; i++)
				{
					if(unvisited.Contains(nearby[i]))
					{
						next = nearby[i];
						available = true;
						break;
					}
				}
				if(!available)
				{
					next = unvisited.ElementAt(rnd.Next() % unvisited.Count);
				}

				int odds = rnd.Next() % Cities.Length;
				if(odds == 0)
				{
					next = unvisited.ElementAt(rnd.Next() % unvisited.Count);
				}
				
				// add that closest city to the route
				route.Add(next);
				unvisited.Remove(next);
			}
			// greedy solution is complete
			
			return new GSCitizen(route);
		}

		private GSCitizen getRandomSolution()
		{
			int i, swap, temp, count = 0;
			string[] results = new string[3];
			int[] perm = new int[Cities.Length];
			ArrayList solution = new ArrayList();
			TSPSolution tspSolution;
			List<int> route = new List<int>();

			do
			{
				for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
					perm[i] = i;
				for (i = 0; i < perm.Length; i++)
				{
					swap = i;
					while (swap == i)
						swap = rnd.Next(0, Cities.Length);
					temp = perm[i];
					perm[i] = perm[swap];
					perm[swap] = temp;
				}
				solution.Clear();
				for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
				{
					solution.Add(Cities[perm[i]]);
					route.Add(perm[i]);
				}
				tspSolution = new TSPSolution(solution);
				count++;
			} while (tspSolution.costOfRoute() == double.PositiveInfinity);                // until a valid route is found

			return new GSCitizen(route);
		}

		public void Shuffle<T>(List<T> list)
		{
			int n = list.Count;
			//Random rnd = new Random();
			while (n > 1)
			{
				int k = (rnd.Next(0, n) % n);
				n--;
				T value = list[k];
				list[k] = list[n];
				list[n] = value;
			}
		}
        #endregion
    }
}
