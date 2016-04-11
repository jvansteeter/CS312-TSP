using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TSP;

namespace WindowsFormsApplication1
{
	class GSCitizen : IComparable<GSCitizen>
	{
		private TSPSolution solution;
		private List<int> route;
		public int[] links;
		private static City[] Cities;
		private static Random rnd;
		private static int citizenCount = 0;

		public static void setCities(City[] Cities)
		{
			GSCitizen.Cities = Cities;
		}

		public static void setRandom(Random rnd)
		{
			GSCitizen.rnd = rnd;
		}

		public GSCitizen(List<int> route)
		{
			GSCitizen.citizenCount++;
			this.route = route;
			this.links = new int[Cities.Length];

			ArrayList solution = new ArrayList();
			for (int i = 0; i < route.Count; i++)
			{
				if (i + 1 == route.Count)
					links[route[i]] = route[0];
				else
					links[route[i]] = route[i + 1];
				solution.Add(Cities[route[i]]);
			}
			this.solution = new TSPSolution(solution);
		}

		/*public GSCitizen(int[] links)
		{
			GSCitizen.citizenCount++;

			this.links = links;
			ArrayList solution = new ArrayList();
			List<int> route = new List<int>();
			
			// construct the TSPSolution from the links
			// first build a route on ints
			for (int i = 0; i < links.Length; i++)
			{
				if (i == 0)
					route.Add(links[0]);
				else
					route.Add(links[route[route.Count - 1]]);
			}
			// using the route, construct the TSPSolution
			for(int i = 0; i < route.Count; i++)
			{
				solution.Add(Cities[route[i]]);
			}
			this.route = route;
			this.solution = new TSPSolution(solution);
		}*/

		public GSCitizen[] reproduce(GSCitizen mate)
		{
			GSCitizen[] children = new GSCitizen[2];
			List<int> route1 = new List<int>();
			List<int> route2 = new List<int>();
			// there is a chance that reproducing could create incomplete cycles, to prevent this keep track of what cities are already linked to
			HashSet<int> missingLinks1 = new HashSet<int>();
			HashSet<int> missingLinks2 = new HashSet<int>();
			for (int i = 0; i < Cities.Length; i++)
			{
				missingLinks1.Add(i);
				missingLinks2.Add(i);
			}

			// breed the two solutions by alternating which parent a child gets a link from
			// do whatever is necessary to prevent multiple cycles
			for(int i = 0; i < GSCitizen.Cities.Length; i++)
			{
				if(i % 2 == 0)
				{
					if(i == 0)
					{
						route1.Add(this.route[0]);
						missingLinks1.Remove(this.route[0]);
						route2.Add(mate.route[0]);
						missingLinks2.Remove(mate.route[0]);
					}
					else 
					{
						if(missingLinks1.Contains(this.links[route1[route1.Count - 1]]))
						{
							int next = this.links[route1[route1.Count - 1]];
							route1.Add(next);
							missingLinks1.Remove(next);
						}
						else
						{
							int next = missingLinks1.ElementAt(rnd.Next() % missingLinks1.Count);
							route1.Add(next);
							missingLinks1.Remove(next);
						}
						if (missingLinks2.Contains(mate.links[route2[route2.Count - 1]]))
						{
							int next = mate.links[route2[route2.Count - 1]];
							route2.Add(next);
							missingLinks2.Remove(next);
						}
						else
						{
							int next = missingLinks2.ElementAt(rnd.Next() % missingLinks2.Count);
							route2.Add(next);
							missingLinks2.Remove(next);
						}
					}
				}
				else
				{
					if (missingLinks1.Contains(mate.links[route1[route1.Count - 1]]))
					{
						int next = mate.links[route1[route1.Count - 1]];
						route1.Add(next);
						missingLinks1.Remove(next);
					}
					else
					{
						int next = missingLinks1.ElementAt(rnd.Next() % missingLinks1.Count);
						route1.Add(next);
						missingLinks1.Remove(next);
					}
					if (missingLinks2.Contains(this.links[route2[route2.Count - 1]]))
					{
						int next = this.links[route2[route2.Count - 1]];
						route2.Add(next);
						missingLinks2.Remove(next);
					}
					else
					{
						int next = missingLinks2.ElementAt(rnd.Next() % missingLinks2.Count);
						route2.Add(next);
						missingLinks2.Remove(next);
					}
				}
			}

			// generate children and give a chance for them to mutate
			GSCitizen child1 = new GSCitizen(route1);
			int chance = rnd.Next() % 100;
			if(chance < 3)
				child1.mutate();
			GSCitizen child2 = new GSCitizen(route2);
			chance = rnd.Next() % 100;
			if(chance < 3)
				child2.mutate();
			children[0] = child1;
			children[1] = child2;
			
			return children;
		}

		public void mutate()
		{
			//Console.WriteLine("Mutating a horrible monster!");
			// Randomly relocate a city in the route
			int randCityIndex = rnd.Next() % this.route.Count;
			int randCity = this.route[randCityIndex];
			this.route.Remove(randCity);
			int randLocation = rnd.Next() % this.route.Count;
			this.route.Insert(randLocation, randCity);
			
			// construct a new TSPSolution from the mutated route
			ArrayList solution = new ArrayList();
			for(int i = 0; i < route.Count; i++)
			{
				solution.Add(Cities[route[i]]);
			}
			TSPSolution tspSolution = new TSPSolution(solution);
			this.solution = tspSolution;
			// create new mutated child and replace this child
			for (int i = 0; i < this.route.Count; i++)
			{
				if (i + 1 == this.route.Count)
					links[this.route[i]] = this.route[0];
				else
					links[this.route[i]] = this.route[i + 1];
			}
		}

		public TSPSolution getSolution()
		{
			return this.solution;
		}

		public double fitness()
		{
			return solution.costOfRoute();
		}

		public bool isValid()
		{
			if (fitness() < Double.PositiveInfinity)
				return true;
			else
				return false;
		}

		public int CompareTo(GSCitizen other)
		{
			if (this.fitness() < other.fitness())
				return -1;
			else
				return 1;
		}
	}
}
