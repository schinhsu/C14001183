import java.util.List;
import java.util.Random;
import java.util.ArrayList;
 
public class Quicksort {
	
	public static void main(String[] args){
		int[] to_sort = new int[10];
		Random ran = new Random();
		for (int i=0;i<10;i++)
			to_sort[i] = ran.nextInt(100)+1;
		System.out.printf("Before sort: ");
		for (int i=0;i<10;i++)
			System.out.printf("%d ",to_sort[i]);
		System.out.printf("\n");
		Sort(to_sort);
		System.out.printf("After sort: ");
		for (int i=0;i<10;i++)
			System.out.printf("%d ",to_sort[i]);
		System.out.printf("\n");
	}
 
    public static void Sort(int[] array){
        List<Integer> list = new ArrayList<Integer>();
        for (int n : array)
            list.add(n);
        list = Sort(list);
        for (int i=0;i<array.length;i++)
            array[i] = list.get(i);
    }
 
    public static List<Integer> Sort(List<Integer> list){
        if (list.size()<2)
            return list;
		
        // middle pivot
        int pivot = list.get(list.size() / 2);
        list.remove(list.size() / 2);
        List<Integer> less = new ArrayList<Integer>();
        List<Integer> greater = new ArrayList<Integer>();
        List<Integer> result = new ArrayList<Integer>();
        for (Integer n : list){
            if (n > pivot) greater.add(n);
            else less.add(n);
        }
        result.addAll(Sort(less));
        result.add(pivot);
        result.addAll(Sort(greater));
        return result;
    }
}