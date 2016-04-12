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
			string[] results = new string[3];
			int solutionCount = 0;
			Stopwatch timer = new Stopwatch();
			PriorityQueue<State> pQueue = new PriorityQueue<State>();
			int maxStoredStates = 0;
			int bssfUpdates = 0;
			int totalStatesCreated = 0;
			int totalStatesPruned = 0;

			timer.Start();
			string[] intialResults = defaultSolveProblem();
			double bestCostSoFar = Int32.Parse(intialResults[COST]);
			Console.WriteLine("BestCostSoFar: " + bestCostSoFar);

			State initialState = initializeState();
			this.reduceState(ref initialState);
			pQueue.insert(initialState, (int)initialState.lowerBound);

			while (!pQueue.isEmpty() && timer.Elapsed.Seconds < 60) // while queue is not empty or has not passed 60 seconds
			{
				State parentState = pQueue.deleteMin();
				if (parentState.lowerBound >= bestCostSoFar) // trim queue if state has worse bound than bssf
				{
					totalStatesPruned++;
					continue;
				}
				for (int i = 0; i < this.Cities.Length; i++)
				{
					if (!parentState.visited.Contains(i))
					{
						// createChild
						totalStatesCreated++;
						State childState = parentState.getCopy();
						this.addCityToRoute(ref childState, i);
						this.reduceState(ref childState);
						if (childState.route.Count == this.Cities.Count()) // if route is complete
						{
							TSPSolution newSolution = new TSPSolution(childState.route);
							double cost = newSolution.costOfRoute();
							if (cost < bestCostSoFar)
							{
								bssfUpdates++;
								bestCostSoFar = cost;
								bssf.Route = childState.route;
							}
						}
						else // route is not complete
						{
							if (childState.lowerBound < bestCostSoFar) // don't insert in queue if state has worse bound than bssf
							{
								pQueue.insert(childState, (int)childState.lowerBound);
								if (pQueue.getCount() > maxStoredStates)
								{
									maxStoredStates = pQueue.getCount();
								}
							}
							else
							{
								totalStatesPruned++;
							}
						}

					}
				}
			}

			timer.Stop();

			results[COST] = costOfBssf().ToString();
			results[TIME] = timer.Elapsed.ToString();
			results[COUNT] = solutionCount.ToString();


			Console.WriteLine("===================================");
			Console.WriteLine("StatesCreated: " + totalStatesCreated);
			Console.WriteLine("StatesPruned: " + totalStatesPruned);
			Console.WriteLine("bssfUpdates: " + bssfUpdates);
			Console.WriteLine("MaxStoredStates: " + maxStoredStates);


			return results;
		}

		private struct State
		{
			public double[,] costMatrix;
			public double lowerBound;
			public ArrayList route;
			public HashSet<int> visited;
			public int lastVisitedIndex;

			public State(double[,] matrix, double lowerBound)
			{
				this.costMatrix = matrix;
				this.lowerBound = lowerBound;
				this.route = new ArrayList();
				this.visited = new HashSet<int>();
				this.lastVisitedIndex = 0;
			}

			public State getCopy()
			{
				State newState = new State();
				newState.costMatrix = (double[,])this.costMatrix.Clone();
				newState.lowerBound = this.lowerBound;
				newState.route = new ArrayList(this.route);
				newState.visited = new HashSet<int>(this.visited);
				newState.lastVisitedIndex = this.lastVisitedIndex;
				return newState;
			}

			public override string ToString()
			{
				StringBuilder sb = new StringBuilder();
				sb.Append(costMatrixToString(this.costMatrix));
				sb.Append("Lowerbound: " + lowerBound + "\n");
				sb.Append("Route: ");
				foreach (var thing in route)
				{
					sb.Append(thing + " ");
				}
				return sb.ToString() + "\n";
			}

			private string costMatrixToString(double[,] costMatrix)
			{
				StringBuilder sb = new StringBuilder();
				for (int i = 0; i < costMatrix.GetLength(0); i++)
				{
					sb.Append("\n");
					for (int j = 0; j < costMatrix.GetLength(1); j++)
					{
						if (costMatrix[i, j] == double.PositiveInfinity)
						{
							sb.Append("~");
						}
						sb.Append(costMatrix[i, j]);
						if (costMatrix[i, j] == double.PositiveInfinity)
						{
							sb.Append("~\t\t");
						}
						else
						{
							sb.Append("\t\t");
						}
					}
				}
				sb.Append("\n");
				return sb.ToString();
			}

		}

		private class PriorityQueue<T>
		{
			private SortedDictionary<int, HashSet<T>> items = new SortedDictionary<int, HashSet<T>>();
			private int count = 0;
			// Total Time Complexity = O(1)
			public void insert(T item, int priority)
			{
				count++;
				if (items.ContainsKey(priority))
				{
					items[priority].Add(item);
				}
				else
				{
					HashSet<T> newSet = new HashSet<T>();
					newSet.Add(item);
					items[priority] = newSet;
				}
			}

			public int getCount()
			{
				return count;
			}

			public T deleteMin()
			{
				count--;
				int key = this.items.Keys.First(); // get key for first set
				T item = this.items[key].First(); // get first item in first set
				this.items[key].Remove(item); // remove item from set
				if (this.items[key].Count == 0) // if set is empty
				{
					this.items.Remove(key); // remove set from dictionary
				}
				return item;
			}

			// Total Time Complexity = O(1)
			public bool isEmpty()
			{
				return (this.items.Count == 0);
			}

			public void printQueue()
			{
				if (this.items.Count == 0)
				{
					Console.WriteLine("Priority Queue is empty!");
					return;
				}
				foreach (int key in this.items.Keys)
				{
					HashSet<T> set = this.items[key];
					Console.WriteLine("Priority: " + key + "\n\t");
					foreach (T item in set)
					{
						Console.Write(item);
					}
				}
			}
		}

		private State initializeState()
		{
			int NUM_CITIES = this.Cities.Length;
			double[,] costMatrix = new double[NUM_CITIES, NUM_CITIES];
			for (int i = 0; i < costMatrix.GetLength(0); i++)
			{
				for (int j = 0; j < costMatrix.GetLength(1); j++)
				{
					if (i == j) // diagonal
					{
						costMatrix[i, j] = double.PositiveInfinity;
					}
					else
					{
						costMatrix[i, j] = this.Cities[i].costToGetTo(this.Cities[j]);
					}
				}
			}

			State newState = new State(costMatrix, 0);
			newState.route.Add(this.Cities[0]);
			newState.visited.Add(0);
			return newState;
		}

		private void addCityToRoute(ref State state, int cityIndex)
		{
			state.lowerBound += state.costMatrix[state.lastVisitedIndex, cityIndex];
			for (int i = 0; i < state.costMatrix.GetLength(0); i++)
			{
				state.costMatrix[i, cityIndex] = double.PositiveInfinity;
			}
			for (int j = 0; j < state.costMatrix.GetLength(1); j++)
			{
				state.costMatrix[state.lastVisitedIndex, j] = double.PositiveInfinity;
			}
			state.costMatrix[cityIndex, state.lastVisitedIndex] = double.PositiveInfinity;
			state.route.Add(this.Cities[cityIndex]);
			state.visited.Add(cityIndex);
			state.lastVisitedIndex = cityIndex;
		}

		private void reduceState(ref State state)
		{
			reduceRows(ref state);
			reduceColumns(ref state);
		}

		private void reduceRows(ref State state)
		{
			for (int i = 0; i < state.costMatrix.GetLength(0); i++)
			{
				double min = state.costMatrix[i, 0]; // first thing in row
				for (int j = 0; j < state.costMatrix.GetLength(1); j++)
				{
					if (state.costMatrix[i, j] < min)
					{
						min = state.costMatrix[i, j];
					}
				}
				if (min > 0 && min != double.PositiveInfinity)
				{
					state.lowerBound += min;
					// substract everything in that row
					for (int j = 0; j < state.costMatrix.GetLength(1); j++)
					{
						state.costMatrix[i, j] = state.costMatrix[i, j] - min;
					}
				}
			}
		}

		private void reduceColumns(ref State state)
		{
			for (int j = 0; j < state.costMatrix.GetLength(0); j++)
			{
				double min = state.costMatrix[0, j]; // first thing in column
				for (int i = 0; i < state.costMatrix.GetLength(1); i++)
				{
					if (state.costMatrix[i, j] < min)
					{
						min = state.costMatrix[i, j];
					}
				}
				if (min > 0 && min != double.PositiveInfinity)
				{
					state.lowerBound += min;
					// substract everything in that column
					for (int i = 0; i < state.costMatrix.GetLength(1); i++)
					{
						state.costMatrix[i, j] = state.costMatrix[i, j] - min;
					}
				}
			}
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
			cityMap = new CityMap(Cities, Cities.Length);

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
					int[] nearest = cityMap.getClosestCities(route[route.Count - 1]);
					int next;
					int count = 0;
					do
					{
						next = nearest[count];
						count++;
					} while (!unvisited.Contains(next));
					unvisited.Remove(next);
					route.Add(next);
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
			int attempts = 0;
			// start time
			Stopwatch timer = new Stopwatch();
			timer.Start();
			// Best Citizen So Far
			GSCitizen bcsf = getRandomSolution();
			solutionsCount++;
			String bestTime = timer.Elapsed.ToString();
			// set global data for Genetic Solution
			int populationSize = 2000;// Cities.Length * 15;
			int groupSize = 5;
			int numNearbyCities = Cities.Length;
            int bestSoFarGen = 0;
			cityMap = new CityMap(Cities, numNearbyCities);
			GSCitizen.cityMap = cityMap;
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
			//while(generationCount < 200)
			while(attempts < populationSize*10)
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
						attempts++;
						// if any survivers are more fit than the best citizen so far, make it the new best citizen so far
						if (survivors[i].fitness() < bcsf.fitness())
						{
							bcsf = survivors[i];
							bestTime = timer.Elapsed.ToString();
							solutionsCount++;
							attempts = 0;
						}
					}
				}
				population = nextGeneration;
				generationCount++;
				//Console.WriteLine("Finished a generation");
			}

			//bssf = getGreedySolution(0).getSolution();
			bssf = bcsf.getSolution();
			double[] genValues = new double[populationSize];
			for(int i = 0; i < populationSize; i++)
			{
				genValues[i] = population[i].fitness();
			}

            //            bestTime = timer.Elapsed.ToString();
            Console.WriteLine("Run time: " + timer.Elapsed.ToString());
            string[] results = new string[3];
            results[COST] = costOfBssf().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
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
