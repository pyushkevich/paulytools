#include "BinaryHeap.h"

int TestBinaryHeapInteractive(int argc, char *argv[])
{
  double weight[] = 
    { 5., 2., 17., 23., 8., 11., 7., 12., 1., 6.};

  BinaryHeap<double> heap(10, weight);

  while(1) {
    string buffer;
    
    cout << " Enter command [i,m,+,-] > " << flush;
    cin >> buffer;

    if(buffer[0] == 'i')
      {
      int iPos = atoi(buffer.c_str()+1);
      heap.InsertElement(iPos);
      }
    else if(buffer[0] == '+')
      {
      int iPos = atoi(buffer.c_str()+1);
      heap.IncreaseElementWeight(iPos, weight[iPos]+1);
      }
    else if(buffer[0] == '-')
      {
      int iPos = atoi(buffer.c_str()+1);
      heap.DecreaseElementWeight(iPos, weight[iPos]-1);
      }
    else if(buffer[0] == 'm')
      {
      cout << "Popped " << heap.PopMinimum() << endl;
      }
    else
      {
      return 0;
      }

    heap.PrintHeap();
  }
}
