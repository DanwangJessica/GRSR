

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;

public class getAInts {

	public static void main(String[] args) throws IOException {
		//args[0] is the position for file 1 which can be empty
		//args[1] is the position for file 2 which can be empty
		//args[2] is the repeat length cutoff
		//output all the start end of the intersection
		//output null if no intersection
		int cutoff=Integer.parseInt(args[2]);
		
		FileReader fr1 = new FileReader(args[0]); 
		BufferedReader br1 = new BufferedReader(fr1);
		String line;
		String[] tokens;
		ArrayList<Block> blks1=new ArrayList<Block>();
		while(true){
			line=br1.readLine();
			if(line==null) break;
			tokens=line.split(" ");
			blks1.add(new Block(Integer.parseInt(tokens[0]),Integer.parseInt(tokens[1])));
		}
		
		FileReader fr2 = new FileReader(args[1]); 
		BufferedReader br2 = new BufferedReader(fr2);
		ArrayList<Block> blks2=new ArrayList<Block>();
		while(true){
			line=br2.readLine();
			if(line==null) break;
			tokens=line.split(" ");
			blks2.add(new Block(Integer.parseInt(tokens[0]),Integer.parseInt(tokens[1])));
		}
		
		if(blks1.size()==0 || blks2.size()==0) return; //if there is empty input file, print nothing
		
		blks1.sort(new CustomComparator());
		//System.out.println("File 1:Before");for(int j=0;j<blks1.size();j++){System.out.println(blks1.get(j).start+" "+blks1.get(j).end);}
		Iterator<Block> it1 = blks1.iterator();
		Block blk1=it1.next();
		Block prev=blk1;
		Block curt;
		while(it1.hasNext()){
			curt=it1.next();
			if(prev.covers(curt))
				it1.remove();
			else
				prev=curt;
		}
		
		//System.out.println("File 1:After");
		//for(int j=0;j<blks1.size();j++){System.out.println(blks1.get(j).start+" "+blks1.get(j).end);}
		
		
		
		blks2.sort(new CustomComparator());
		//System.out.println("File 2:Before");for(int j=0;j<blks2.size();j++){System.out.println(blks2.get(j).start+" "+blks2.get(j).end);}
		Iterator<Block> it2 = blks2.iterator();
		Block blk2=it2.next();
		prev=blk2;
		while(it2.hasNext()){
			curt=it2.next();
			if(prev.covers(curt))
				it2.remove();
			else
				prev=curt;
		}
		//System.out.println("File 2:After");
		//for(int j=0;j<blks2.size();j++){System.out.println(blks2.get(j).start+" "+blks2.get(j).end);}
		
		
		
		int sp,ep,s,e; //sp,ep (s,p) is the start, end for the blocks in blks2 (blks1); 
		int li,si,ei; //length, start, end of the intersection segment
		
		//System.out.println("Intersection:");
		for(int i=0;i<blks2.size();i++){
			sp=blks2.get(i).start;
			ep=blks2.get(i).end; //block prime in blks2;
			
			for(int j=0;j<blks1.size();j++){
				s=blks1.get(j).start;
				e=blks1.get(j).end;
			
				//The remaining is for the two cases which no intersection between this two blocks
				if(e<=sp) //block prime is after block
					continue;
				if(ep<=s)//block prime is before block
					break;
				
				//The remaining is for cases which intersection must exist
				si=Math.max(s,sp);
				ei=Math.min(e,ep);
				li=ei-si+1;
				if(li>=cutoff)
					System.out.println(si+" "+ei);
		
			}
		}
		
		br1.close();
		br2.close();

	}

}

class Block{
	int start;
	int end;
	public Block(int s,int e){
		this.start=s;
		this.end=e;
	}
	public boolean covers(Block b){
		if(start<=b.start && end>=b.end){
			return true;
		}
		else
			return false;
	}
}

class CustomComparator implements Comparator<Block> {
    public int compare(Block o1, Block o2) {
        return o1.start - o2.start;
    }
}


