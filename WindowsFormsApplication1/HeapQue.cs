using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TSP;

namespace WindowsFormsApplication1
{
    class HeapQue
    {
        private double[] dist;
        private int[] heap;
        private int[] pointer;
        private int[] prev;
        private int endIndex;

        public HeapQue(List<City> points)
        {
            // Make Que
            // Size O(n)
            // Time O(n)
            dist = new double[points.Count];
            pointer = new int[points.Count];
            heap = new int[points.Count];
            prev = new int[points.Count];
            endIndex = points.Count - 1;

            for (int i = 0; i < points.Count; i++)
            {
                dist[i] = Double.MaxValue;
                heap[i] = i;
                pointer[i] = i;
                prev[i] = -1;
            }
        }

        public int deleteMin()
        {
            // Space O(1)
            // Time O(logn)
            int minIndex = heap[0];
            pointer[heap[0]] = -1;
            heap[0] = heap[endIndex];
            pointer[heap[0]] = 0;
            heap[endIndex] = -1;
            endIndex--;

            int currentIndex = 0;
            while (currentIndex < endIndex)
            {
                int child1 = (currentIndex + 1) * 2;
                int child2 = ((currentIndex + 1) * 2) - 1;

                if (child1 >= heap.Length || heap[child1] == -1)
                    break;

                if ((child1 < heap.Length && child2 >= heap.Length) || (child1 < heap.Length && heap[child2] == -1))
                {
                    if (dist[heap[currentIndex]] > dist[heap[child1]])
                    {
                        int temp = heap[currentIndex];
                        heap[currentIndex] = heap[child1];
                        pointer[heap[currentIndex]] = currentIndex;
                        heap[child1] = temp;
                        pointer[heap[child1]] = child1;
                    }
                    break;
                }

                if (dist[heap[currentIndex]] > dist[heap[child1]] || dist[heap[currentIndex]] > dist[heap[child2]])
                {
                    if (dist[heap[child1]] < dist[heap[child2]])
                    {
                        int temp = heap[currentIndex];
                        heap[currentIndex] = heap[child1];
                        pointer[heap[currentIndex]] = currentIndex;
                        heap[child1] = temp;
                        pointer[heap[child1]] = child1;
                        currentIndex = child1;
                    }
                    else
                    {
                        int temp = heap[currentIndex];
                        heap[currentIndex] = heap[child2];
                        pointer[heap[currentIndex]] = currentIndex;
                        heap[child2] = temp;
                        pointer[heap[child2]] = child2;
                        currentIndex = child2;
                    }
                }
                else
                    break;
            }
            return minIndex;
        }

        public void decreaseKey(int index, double value)
        {
            // Space O(1)
            // Time O(1)
            int inHeap = pointer[index];
            dist[index] = value;

            int currentIndex = inHeap;
            while (currentIndex > -1)
            {
                int parent = (currentIndex - 1) / 2;
                if (dist[heap[currentIndex]] < dist[heap[parent]])
                {
                    int temp = heap[currentIndex];
                    heap[currentIndex] = heap[parent];
                    pointer[heap[currentIndex]] = currentIndex;
                    heap[parent] = temp;
                    pointer[heap[parent]] = parent;
                    currentIndex = parent;
                }
                else
                {
                    break;
                }
            }
        }

        public bool isEmpty()
        {
            // Space O(1)
            // Time O(1)
            if (endIndex == -1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public double getDistance(int index)
        {
            // Space O(1)
            // Time O(1)
            return dist[index];
        }

        public void setPrev(int index, int value)
        {
            // Space O(1)
            // Time O(1)
            prev[index] = value;
        }

        public int getPrev(int index)
        {
            // Space O(1)
            // Time O(1)
            return prev[index];
        }
    }
}
