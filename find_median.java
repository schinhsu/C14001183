import java.util.List;
import java.util.Random;
import java.util.ArrayList;
 
public class find_median {
	public static void main(String[] args){
		int[] to_sort = new int[9];
		Random ran = new Random();
		for (int i=0;i<9;i++)
			to_sort[i] = ran.nextInt(100)+1;
		System.out.printf("Before sort: ");
		for (int i=0;i<9;i++)
			System.out.printf("%d ",to_sort[i]);
		System.out.printf("\n");
		Partition(to_sort);
		System.out.printf("Find median: %d\n",to_sort[0]);
	}
 
    public static void Partition(int[] array){
        List<Integer> list = new ArrayList<Integer>();
        for (int n : array)
            list.add(n);
        list = Partition(list,(list.size()+1)/2);
        array[0] = list.get(0);
    }
 
    public static List<Integer> Partition(List<Integer> list,int index){
        if (list.size()<2)
            return list;
		
        // middle pivot
        int pivot = list.get(list.size() / 2);
		int l1 = 0;
        list.remove(list.size() / 2);
        List<Integer> less = new ArrayList<Integer>();
        List<Integer> greater = new ArrayList<Integer>();
        List<Integer> result = new ArrayList<Integer>();
        for (Integer n : list){
            if (n > pivot){
				greater.add(n);
			}
            else{
				less.add(n);
				l1++;
			}
        }
		if (index == l1+1){
			result.add(pivot);
			return result;
		}
		else if (index < l1+1){
			result.addAll(less);
		}
		else{
			result.addAll(greater);
			index = index-l1-1;
		}
		return Partition(result,index);
    }
}
