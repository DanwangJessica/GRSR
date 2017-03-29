
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

public class deleTransp {
	public static void main(String[] args) throws IOException {
		//args[0] is the number of blocks
		//args[1~n-1] is the source and destination scaffold
		
		ArrayList<ArrayList<Integer>> genomePair=new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<2;i++) genomePair.add(new ArrayList<Integer>());
		
		int noBlk=Integer.parseInt(args[0]);
		String g0="";
		String g1="";
		for(int i=1;i<=noBlk;i++){
			g0+=args[i]+" ";
			g1+=args[i+noBlk+1]+" ";
		}
		g0+="$";
		g1+="$";
		
		genomePair=inOrder(g0,g1,"BefMerge");
		genomePair=merge2Blks(genomePair,"AftMerge");
		
		g0=genomePair.get(0).toString().replace("[","").replace("]","").replace(",","")+" $";
		g1=genomePair.get(1).toString().replace("[","").replace("]","").replace(",","")+" $";
		genomePair=inOrder(g0,g1,"BefFilter");
		
		genomePair=filterTBI(genomePair);
		
		g0=genomePair.get(0).toString().replace("[","").replace("]","").replace(",","")+" $";
		g1=genomePair.get(1).toString().replace("[","").replace("]","").replace(",","")+" $";
		genomePair=inOrder(g0,g1,"AftFilter");
		
		PrintWriter writer=new PrintWriter("mgr_macro2"+".txt","UTF-8");
		g0=genomePair.get(0).toString().replace("[","").replace("]","").replace(",","")+" $";
		g1=genomePair.get(1).toString().replace("[","").replace("]","").replace(",","")+" $";
		writer.println(">genome1");
		writer.println(g0);
		writer.println(">genome2");
		writer.println(g1);
		writer.close();
	}
	//This function is for make the genomes in the correct order, s and d should ended with $
	public static ArrayList<ArrayList<Integer>> inOrder(String s,String d,String filename) throws FileNotFoundException, UnsupportedEncodingException{
		ArrayList<ArrayList<Integer>> orderedGenomes=new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<2;i++) orderedGenomes.add(new ArrayList<Integer>());
		PrintWriter writer=new PrintWriter(filename+".txt","UTF-8");
		
		String line=s; //get the first genome's scaffold
		line=line.replace(" $","");
		String[] sGenome=line.split(" ");
		String line1=d; //get the second genome's scaffold
		line1=line1.replace(" $","");
		String[] dGenome=line1.split(" ");
		
		int blkNo=sGenome.length;
		String oldBlk;
		int newBlk; //store the new blkID for destination
		int id;
		int newsBlk; //store the new blkID for source
		for(int h=0;h<blkNo;h++){
			newsBlk=h+1;
			orderedGenomes.get(0).add(newsBlk);
			writer.println(sGenome[h]+"="+newsBlk); //old_blockID:new_blockID
			oldBlk=dGenome[h];
			id=Arrays.asList(sGenome).indexOf(oldBlk); // current block is the same as those in source genome
			if(id>=0) 
				newBlk=id+1;
			else{ //destination and source, this blk is in reverse direction
				id=Arrays.asList(sGenome).indexOf("-"+oldBlk); //d=3 s=-3
				if(id>=0) 
					newBlk=(-1)*(id+1);
				else{ //d=-3 s=3
					id=Arrays.asList(sGenome).indexOf(oldBlk.replace("-",""));
					newBlk=(-1)*(id+1);
				}
			}
			orderedGenomes.get(1).add(newBlk);
		}
		
		writer.close();
		return orderedGenomes;
	}
	//This function is for merging pairwised genomes
	public static ArrayList<ArrayList<Integer>> merge2Blks(ArrayList<ArrayList<Integer>> genomePair,String filename) throws IOException{
		final int NO_OF_GENOMES=2;
		int blkNo=genomePair.get(0).size();
		ArrayList<ArrayList<Integer>> genomes=new ArrayList<ArrayList<Integer>>();
		for(int c=0;c<NO_OF_GENOMES;c++) genomes.add(new ArrayList<Integer>());
		int onePos,monePos;
		int c1;
		for(int j=0;j<NO_OF_GENOMES;j++){
			onePos=genomePair.get(j).indexOf(1);
			if(onePos==0){
				for(int k=0;k<blkNo;k++) genomes.get(j).add(genomePair.get(j).get(k));
			}
			else if(onePos>0){
				c1=0;
				for(int k=0;k<blkNo;k++){
					if(onePos+k<blkNo)
						genomes.get(j).add(genomePair.get(j).get(onePos+k));
					else{
						genomes.get(j).add(genomePair.get(j).get(c1));
						c1++;
					}
				}
			}
			else{ //means find the -1 in strain j
				monePos=genomePair.get(j).indexOf(-1);
				c1=0;
				for(int k=0;k<blkNo;k++){
					if(monePos-k>=0)
						genomes.get(j).add((-1)*genomePair.get(j).get(monePos-k));
					else{
						genomes.get(j).add((-1)*genomePair.get(j).get(blkNo-1-c1));
						c1++;
					}
				}
			}
		}
		
		Iterator<Integer> it = genomes.get(0).iterator();
		int cblk;
		int[] p=new int[NO_OF_GENOMES];
		Integer g0nextB;
		Boolean allEq=true;
		int[] count=new int[NO_OF_GENOMES];
		
		while(it.hasNext()){
			cblk=it.next();
			//search the position of current block in other 24 genomes
			for(int j1=0;j1<NO_OF_GENOMES;j1++){ 
				p[j1]=genomes.get(j1).indexOf(cblk); //if find the current block itself, the position is >=0
				if(p[j1]<0) p[j1]=(-1)*genomes.get(j1).indexOf(-cblk); //if find the current block's negative value, the position is <0.
				
			}
			
			//check whether the next block are all equal
			if(it.hasNext()){
				g0nextB=genomes.get(0).get(p[0]+1);
			}
			else
				break;
			
			for(int k=1;k<NO_OF_GENOMES;k++){
				if(p[k]>=0 && p[k]<genomes.get(k).size()-1 && !genomes.get(k).get(p[k]+1).equals(g0nextB)){
					allEq=false;
					break;
				}
				if(p[k]==genomes.get(k).size()-1 && !genomes.get(k).get(0).equals(g0nextB)){
					
					allEq=false;
					break;
				}
				if(p[k]<0 && !genomes.get(k).get(-p[k]-1).equals((-1)*g0nextB)){
					
					allEq=false;
					break;
				}
			}
			
			for(int i1=0;i1<NO_OF_GENOMES;i1++) count[i1]=0; //count the number of removed elements
			while(allEq){
				//remove the all Equal value
				it.next();
				it.remove();
				for(int h=1;h<NO_OF_GENOMES;h++){
					if (p[h]>=0){
						if(p[h]<genomes.get(h).size()-1)
							genomes.get(h).remove(p[h]+1);
						
						else if(p[h]==genomes.get(h).size()-1)
							genomes.get(h).remove(0);
					}	
					
					if(p[h]<0){
						genomes.get(h).remove(-p[h]-count[h]-1);
						count[h]++;
					}
				}

				if(!it.hasNext()) 
					break;
				//check whether the next block are equal in the 25 genomes
				g0nextB=genomes.get(0).get(p[0]+1); 
				
				//System.out.println(g0nextB);
				for(int k=1;k<NO_OF_GENOMES;k++){
					if(p[k]>=0 && p[k]<genomes.get(k).size()-1 && !genomes.get(k).get(p[k]+1).equals(g0nextB)){
						allEq=false;
						break;
					}
					if(p[k]==genomes.get(k).size()-1 && !genomes.get(k).get(0).equals(g0nextB)){
						allEq=false;
						break;
					}
					if(p[k]<0 && !genomes.get(k).get(-p[k]-count[k]-1).equals((-1)*g0nextB)){
						allEq=false;
						break;
					}
				}
			}
			allEq=true;
		}
		PrintWriter writer=new PrintWriter(filename+".txt","UTF-8");
		int end=0;
		for(int m=0;m<genomes.get(0).size();m++){
			if(m!=genomes.get(0).size()-1){
				end=genomes.get(0).get(m+1)-1;
				writer.println(genomes.get(0).get(m)+":"+end+"="+genomes.get(0).get(m));
			}
			else{
				end=genomePair.get(0).get(blkNo-1);
				writer.println(genomes.get(0).get(m)+":"+end+"="+genomes.get(0).get(m));
			}
		}
		writer.close();
		
		return genomes;
	}
	//This function is for filtering out transposition and block interchange
	public static ArrayList<ArrayList<Integer>> filterTBI(ArrayList<ArrayList<Integer>> genomePair) throws IOException{
		
	    
		int blkNo=genomePair.get(0).size();
		ArrayList<ArrayList<Integer>> lists=new ArrayList<ArrayList<Integer>>();
		for(int c=0;c<4;c++) lists.add(new ArrayList<Integer>()); //0,2 for deletion; 1,3 for insertion
		
		for(int i=0;i<=blkNo-2;i++){
			if(genomePair.get(1).get(i).equals(genomePair.get(1).get(i+1)-2))
				lists.get(0).add(genomePair.get(1).get(i)+1);
		}
		for(int i=0;i<=blkNo-3;i++){
			if(genomePair.get(1).get(i).equals(genomePair.get(1).get(i+2)-1))
				lists.get(1).add(genomePair.get(1).get(i+1));
			if( !genomePair.get(1).get(i+1).equals(genomePair.get(1).get(i)+1) && genomePair.get(1).get(i+2).equals(genomePair.get(1).get(i)+2)){
				lists.get(2).add(genomePair.get(1).get(i)+1);
				lists.get(3).add(genomePair.get(1).get(i+1));
			}
		}
		
		//special case for the the last block
		if(blkNo>=2 && genomePair.get(1).get(0).equals(genomePair.get(0).get(0))){
			//transposition on last block
			if( genomePair.get(1).get(blkNo-1).equals(genomePair.get(0).get(blkNo-2))){
				lists.get(0).add(genomePair.get(0).get(blkNo-1));
			}
			
			if( genomePair.get(1).get(blkNo-2).equals(genomePair.get(0).get(blkNo-1)) ){
				lists.get(1).add(genomePair.get(1).get(blkNo-1));
			}
			
			//block interchange
			if( !genomePair.get(1).get(blkNo-1).equals(genomePair.get(0).get(blkNo-1)) && genomePair.get(1).get(blkNo-2).equals(genomePair.get(0).get(blkNo-2))){
				lists.get(2).add(genomePair.get(0).get(blkNo-1));
				lists.get(3).add(genomePair.get(1).get(blkNo-1));
			}
		}
		//delete blocks for transposition and block interchange;
		int tNo=lists.get(1).size();
		
		ArrayList<Integer> tp=new ArrayList<Integer>(); //Blocks who has been transposed
		ArrayList<Integer> Itp=new ArrayList<Integer>(); //Blocks who has been reversely transposed
		ArrayList<Integer> bi=new ArrayList<Integer>(); //Blocks who has been block interchanged
		ArrayList<Integer> Ibi=new ArrayList<Integer>(); //Blocks who has been reversely block interchanged
		ArrayList<Integer> Hibi=new ArrayList<Integer>(); //Blocks who has been half inverted block interchanged
		
		for(int i=0;i<tNo;i++){
			if(lists.get(0).contains(lists.get(1).get(i))){
				tp.add(new Integer(Math.abs(lists.get(1).get(i))));
				genomePair.get(1).remove(lists.get(1).get(i));
				genomePair.get(0).remove(new Integer(Math.abs(lists.get(1).get(i))));
			}
			else if(lists.get(0).contains((-1)*lists.get(1).get(i))){
				Itp.add(new Integer(Math.abs(lists.get(1).get(i))));
				genomePair.get(1).remove(lists.get(1).get(i));
				genomePair.get(0).remove(new Integer(Math.abs(lists.get(1).get(i))));
			}
		}
		

		Iterator<Integer> it2 = lists.get(2).iterator();
		Iterator<Integer> it3 = lists.get(3).iterator();
		Integer c2blk;
		Integer c3blk;
		int id;
		while(it2.hasNext() && it3.hasNext()){
			c2blk=it2.next();
			c3blk=it3.next();
			if(Math.abs(c2blk)==Math.abs(c3blk)){
				it2.remove();
				it3.remove();
				continue;
			}
			
			boolean inv=false;
			id=lists.get(2).indexOf(c3blk);
			if(id<0){
				id=lists.get(2).indexOf((-1)*c3blk); 
				if(id>0) inv=true;
			}
			
			if(id>0){ //found c3blk in list2
				if(lists.get(3).get(id).equals(c2blk)){
					if(inv){
						Hibi.add(new Integer(Math.abs(lists.get(3).get(id))));
						Hibi.add(new Integer(Math.abs(c3blk)));
					}
					else{
						bi.add(new Integer(Math.abs(lists.get(3).get(id))));
						bi.add(new Integer(Math.abs(c3blk)));
					}
					genomePair.get(1).remove(lists.get(3).get(id));
					genomePair.get(1).remove(c3blk);
					genomePair.get(0).remove(new Integer(Math.abs(lists.get(3).get(id))));
					genomePair.get(0).remove(new Integer(Math.abs(c3blk)));
					it2.remove();
					it3.remove();
				}
				else if(lists.get(3).get(id).equals((-1)*c2blk)){
					if(inv){
						Ibi.add(new Integer(Math.abs(lists.get(3).get(id))));
						Ibi.add(new Integer(Math.abs(c3blk)));
					}
					else{
						Hibi.add(new Integer(Math.abs(lists.get(3).get(id))));
						Hibi.add(new Integer(Math.abs(c3blk)));
					}
					genomePair.get(1).remove(lists.get(3).get(id));
					genomePair.get(1).remove(c3blk);
					genomePair.get(0).remove(new Integer(Math.abs(lists.get(3).get(id))));
					genomePair.get(0).remove(new Integer(Math.abs(c3blk)));
					it2.remove();
					it3.remove();
				}
				else
					continue;
			}
			else //not found c3blk in list3
				continue;
		}
		
		//To deal with adjacent blocks which are involved transposition in tp, e.g.9 11 10 12, Block 10 and 11 both involved transposition
		if(tp.size()>1){
			Iterator<Integer> it4 = tp.iterator();
			Integer cur4;
			Integer prev4;
			prev4=it4.next();
			while(it4.hasNext()){
				cur4=it4.next();
				if(Math.abs(cur4.intValue()-prev4.intValue())==1){ //remove 11, just consider block 10. Removing 11 will not change the breakpoint of this transpostion
					it4.remove();
				}
				else
					prev4=cur4;
				
			}
		}
		
		System.out.println(tp.size()+Itp.size()+bi.size()/2+Ibi.size()/2+Hibi.size()/2);
		System.out.println(tp.size()+" "+Itp.size()+" "+bi.size()/2+" "+Ibi.size()/2+" "+Hibi.size()/2);
		if(tp.size()!=0)
			System.out.println("T "+tp.toString().replace("[","").replace("]","").replace(",",""));
		if(Itp.size()!=0)
			System.out.println("IT "+Itp.toString().replace("[","").replace("]","").replace(",",""));
		if(bi.size()!=0)
			System.out.println("B "+bi.toString().replace("[","").replace("]","").replace(",",""));
		if(Ibi.size()!=0)
			System.out.println("IB "+Ibi.toString().replace("[","").replace("]","").replace(",",""));
		if(Hibi.size()!=0)
			System.out.println("HIB "+Hibi.toString().replace("[","").replace("]","").replace(",",""));
		
	
		return genomePair;
	}
	
	
	
}
