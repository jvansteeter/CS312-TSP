using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TSP;

namespace WindowsFormsApplication1
{
	class SpawningPool
	{
		private PriorityQueue<GSCitizen> que;

		public SpawningPool()
		{
			que = new PriorityQueue<GSCitizen>();
		}

		public void Add(GSCitizen citizen)
		{
			que.Enqueue(citizen);
		}

		public GSCitizen[] spawn()
		{
			// determine the most fit citizens in the spawning pool and breed them
			GSCitizen[] survivors = new GSCitizen[que.Count()];
			GSCitizen parent1 = que.Dequeue();
			GSCitizen parent2 = que.Dequeue();

			GSCitizen[] children = parent1.reproduce(parent2);

			// detemine array of survivors to return
			int survivorCount = 0;
			survivors[0] = parent1;
			survivorCount++;
			survivors[1] = parent2;
			survivorCount++;
			if(children[0].isValid())
			{
				survivors[survivorCount] = children[0];
				survivorCount++;
			}
			if(children[1].isValid())
			{
				survivors[survivorCount] = children[1];
				survivorCount++;
			}
			while(survivors.Length > survivorCount)
			{
				survivors[survivorCount] = que.Dequeue();
				survivorCount++;
			}

			return survivors;
		}
	}
}
