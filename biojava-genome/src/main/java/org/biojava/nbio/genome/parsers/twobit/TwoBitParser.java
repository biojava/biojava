/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.genome.parsers.twobit;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.HashMap;

/**
 * downloaded from http://storage.bioinf.fbb.msu.ru/~roman/TwoBitParser.java
 *
 * Class is a parser of UCSC Genome Browser file format .2bit used to store
 * nucleotide sequence information. This class extends InputStream and can
 * be used as it after choosing one of names of containing sequences. This
 * parser can be used to do some work like UCSC tool named twoBitToFa. For
 * it just run this class with input file path as single parameter and set
 * stdout stream into output file. If you have any problems or ideas don't
 * hesitate to contact me through email: rsutormin[at]gmail.com.
 * @author Roman Sutormin
 */
public class TwoBitParser extends InputStream {

	private static final Logger logger = LoggerFactory.getLogger(TwoBitParser.class);

	public int DEFAULT_BUFFER_SIZE = 10000;
	//
	private RandomAccessFile raf;
	private File f;
	private boolean reverse = false;
	private String[] seq_names;
	private HashMap<String,Long> seq2pos = new HashMap<String,Long>();
	private String cur_seq_name;
	private long[][] cur_nn_blocks;
	private long[][] cur_mask_blocks;
	private long cur_seq_pos;
	private long cur_dna_size;
	private int cur_nn_block_num;
	private int cur_mask_block_num;
	private int[] cur_bits;
	private byte[] buffer;
	private long buffer_size;
	private long buffer_pos;
	private long start_file_pos;
	private long file_pos;
	//
	private static final char[] bit_chars = {
			'T','C','A','G'
	};

	public TwoBitParser(File f) throws Exception {
		this.f = f;
		raf = new RandomAccessFile(f,"r");
		long sign = readFourBytes();
		if(sign==0x1A412743) {
			logger.debug("2bit: Normal number architecture");
		}
		else if(sign==0x4327411A) {
			reverse = true;
			logger.debug("2bit: Reverse number architecture");
		}
		else throw new Exception("Wrong start signature in 2BIT format");
		readFourBytes();
		int seq_qnt = (int)readFourBytes();
		readFourBytes();
		seq_names = new String[seq_qnt];
		for(int i=0;i<seq_qnt;i++) {
			int name_len = raf.read();
			char[] chars = new char[name_len];
			for(int j=0;j<name_len;j++) chars[j] = (char)raf.read();
			seq_names[i] = new String(chars);
			long pos = readFourBytes();
			seq2pos.put(seq_names[i],pos);
			logger.debug("2bit: Sequence name=[{}], pos={}", seq_names[i], pos);
		}
	}
	private long readFourBytes() throws Exception {
		long ret = 0;
		if(!reverse) {
			ret = raf.read();
			ret += raf.read()*0x100;
			ret += raf.read()*0x10000;
			ret += raf.read()*0x1000000;
		}
		else {
			ret = raf.read()*0x1000000;
			ret += raf.read()*0x10000;
			ret += raf.read()*0x100;
			ret += raf.read();
		}
		return ret;
	}
	public String[] getSequenceNames() {
		String[] ret = new String[seq_names.length];
		System.arraycopy(seq_names,0,ret,0,seq_names.length);
		return ret;
	}
	/**
	 * Method open nucleotide stream for sequence with given name.
	 * @param seq_name name of sequence (one of returned by getSequenceNames()).
	 * @throws Exception
	 */
	public void setCurrentSequence(String seq_name) throws Exception {
		if(cur_seq_name!=null) {
			throw new Exception("Sequence ["+cur_seq_name+"] was not closed");
		}
		if(seq2pos.get(seq_name)==null) {
			throw new Exception("Sequence ["+seq_name+"] was not found in 2bit file");
		}
		cur_seq_name = seq_name;
		long pos = seq2pos.get(seq_name);
		raf.seek(pos);
		long dna_size = readFourBytes();
		logger.debug("2bit: Sequence name=[{}], dna_size={}", cur_seq_name, dna_size);
		cur_dna_size = dna_size;
		int nn_block_qnt = (int)readFourBytes();
		cur_nn_blocks = new long[nn_block_qnt][2];
		for(int i=0;i<nn_block_qnt;i++) {
			cur_nn_blocks[i][0] = readFourBytes();
		}
		for(int i=0;i<nn_block_qnt;i++) {
			cur_nn_blocks[i][1] = readFourBytes();
		}

		for(int i=0;i<nn_block_qnt;i++) {
			logger.debug("NN-block: [{},{}] ", cur_nn_blocks[i][0], cur_nn_blocks[i][1]);
		}

		int mask_block_qnt = (int)readFourBytes();
		cur_mask_blocks = new long[mask_block_qnt][2];
		for(int i=0;i<mask_block_qnt;i++) {
			cur_mask_blocks[i][0] = readFourBytes();
		}
		for(int i=0;i<mask_block_qnt;i++) {
			cur_mask_blocks[i][1] = readFourBytes();
		}

		for(int i=0;i<mask_block_qnt;i++) {
			logger.debug("[{},{}] ", cur_mask_blocks[i][0], cur_mask_blocks[i][1]);
		}

		readFourBytes();
		start_file_pos = raf.getFilePointer();
		reset();
	}
	/**
	 * Method resets current position to the begining of sequence stream.
	 */
	@Override
	public synchronized void reset() throws IOException {
		cur_seq_pos = 0;
		cur_nn_block_num = (cur_nn_blocks.length>0)?0:-1;
		cur_mask_block_num = (cur_mask_blocks.length>0)?0:-1;
		cur_bits = new int[4];
		file_pos = start_file_pos;
		buffer_size = 0;
		buffer_pos = -1;
	}
	/**
	 * @return number (starting from 0) of next readable nucleotide in sequence stream.
	 */
	public long getCurrentSequencePosition() {
		if(cur_seq_name==null) throw new RuntimeException("Sequence is not set");
		return cur_seq_pos;
	}
	public void setCurrentSequencePosition(long pos) throws IOException {
		if(cur_seq_name==null) throw new RuntimeException("Sequence is not set");
		if(pos>cur_dna_size) throw new RuntimeException(
				"Postion is too high (more than "+cur_dna_size+")");
		if(cur_seq_pos>pos) {
			reset();
		}
		skip(pos-cur_seq_pos);
	}
	private void loadBits() throws IOException {
		if((buffer==null)||(buffer_pos<0)||(file_pos<buffer_pos)||
				(file_pos>=buffer_pos+buffer_size)) {
			if((buffer==null)||(buffer.length!=DEFAULT_BUFFER_SIZE)) {
				buffer = new byte[DEFAULT_BUFFER_SIZE];
			}
			buffer_pos = file_pos;
			buffer_size = raf.read(buffer);
		}
		int cur_byte = buffer[(int)(file_pos-buffer_pos)]& 0xff;
		for(int i=0;i<4;i++) {
			cur_bits[3-i] = cur_byte%4;
			cur_byte /= 4;
		}
	}
	/**
	 * Method reads 1 nucleotide from sequence stream. You should set current sequence
	 * before use it.
	 */
	@Override
	public int read() throws IOException {
		if(cur_seq_name==null) throw new IOException("Sequence is not set");
		if(cur_seq_pos==cur_dna_size) {
			logger.debug("End of sequence (file position:{})", raf.getFilePointer());
			return -1;
		}
		int bit_num = (int)cur_seq_pos%4;
		if(bit_num==0) {
			loadBits();
		}
		else if(bit_num==3) {
			file_pos++;
		}
		char ret = 'N';
		if((cur_nn_block_num>=0)&&
				(cur_nn_blocks[cur_nn_block_num][0]<=cur_seq_pos)) {
			if(cur_bits[bit_num]!=0) {
				throw new IOException("Wrong data in NN-block ("+cur_bits[bit_num]+") "+
						"at position "+cur_seq_pos);
			}
			if(cur_nn_blocks[cur_nn_block_num][0]+cur_nn_blocks[cur_nn_block_num][1]==cur_seq_pos+1) {
				cur_nn_block_num++;
				if(cur_nn_block_num>=cur_nn_blocks.length) {
					cur_nn_block_num = -1;
				}
			}
			ret = 'N';
		}
		else {
			ret = bit_chars[cur_bits[bit_num]];
		}
		if((cur_mask_block_num>=0)&&
				(cur_mask_blocks[cur_mask_block_num][0]<=cur_seq_pos)) {
			ret = Character.toLowerCase(ret);
			if(cur_mask_blocks[cur_mask_block_num][0]+cur_mask_blocks[cur_mask_block_num][1]==cur_seq_pos+1) {
				cur_mask_block_num++;
				if(cur_mask_block_num>=cur_mask_blocks.length) {
					cur_mask_block_num = -1;
				}
			}
		}
		cur_seq_pos++;
		return ret;
	}
	/**
	 * Method skips n nucleotides in sequence stream. You should set current sequence
	 * before use it.
	 */
	@Override
	public synchronized long skip(long n) throws IOException {
		if(cur_seq_name==null) throw new IOException("Sequence is not set");
		if(n<4) {
			int ret = 0;
			while((ret<n)&&(read()>=0)) ret++;
			return ret;
		}
		if(n>cur_dna_size-cur_seq_pos) {
			n = cur_dna_size-cur_seq_pos;
		}
		cur_seq_pos += n;
		file_pos = start_file_pos+(cur_seq_pos/4);
		raf.seek(file_pos);
		if((cur_seq_pos%4)!=0) {
			loadBits();
		}
		while((cur_nn_block_num>=0)&&
				(cur_nn_blocks[cur_nn_block_num][0]+cur_nn_blocks[cur_nn_block_num][1]<=cur_seq_pos)) {
			cur_nn_block_num++;
			if(cur_nn_block_num>=cur_nn_blocks.length) cur_nn_block_num = -1;
		}
		while((cur_mask_block_num>=0)&&
				(cur_mask_blocks[cur_mask_block_num][0]+cur_mask_blocks[cur_mask_block_num][1]<=cur_seq_pos)) {
			cur_mask_block_num++;
			if(cur_mask_block_num>=cur_mask_blocks.length) cur_mask_block_num = -1;
		}
		return n;
	}
	/**
	 * Method closes current sequence and it's necessary to invoke it before setting
	 * new current sequence.
	 */
	@Override
	public void close() throws IOException {
		cur_seq_name = null;
		cur_nn_blocks = null;
		cur_mask_blocks = null;
		cur_seq_pos = -1;
		cur_dna_size = -1;
		cur_nn_block_num = -1;
		cur_mask_block_num = -1;
		cur_bits = null;
		buffer_size = 0;
		buffer_pos = -1;
		file_pos = -1;
		start_file_pos = -1;
	}
	@Override
	public int available() throws IOException {
		if(cur_seq_name==null) throw new IOException("Sequence is not set");
		return (int)(cur_dna_size-cur_seq_pos);
	}
	/**
	 * Method closes random access file descriptor. You can't use any reading methods
	 * after it.
	 * @throws Exception
	 */
	public void closeParser() throws Exception {
		raf.close();
	}
	public File getFile() {
		return f;
	}
	public String loadFragment(long seq_pos,int len) throws IOException {
		if(cur_seq_name==null) throw new IOException("Sequence is not set");
		setCurrentSequencePosition(seq_pos);
		char[] ret = new char[len];
		int i = 0;
		for(;i<len;i++) {
			int ch = read();
			if(ch<0) break;
			ret[i] = (char)ch;
		}
		return new String(ret,0,i);
	}
	public void printFastaSequence() throws IOException {
		if(cur_seq_name==null) throw new RuntimeException("Sequence is not set");
		printFastaSequence(cur_dna_size-cur_seq_pos);
	}
	public void printFastaSequence(long len) throws IOException {
		if(cur_seq_name==null) throw new RuntimeException("Sequence is not set");
		logger.info(">{} pos={}, len={}", cur_seq_name, cur_seq_pos, len);
		char[] line = new char[60];
		boolean end = false;
		long qnt_all = 0;
		while(!end) {
			int qnt = 0;
			for(;(qnt<line.length)&&(qnt_all<len);qnt++,qnt_all++) {
				int ch = read();
				if(ch<0) {
					end = true;
					break;
				}
				line[qnt] = (char)ch;
			}
			if(qnt>0) {
				logger.info(new String(line,0,qnt));
			}
			if(qnt_all>=len) end = true;
		}
	}
	public static void main(String[] args) throws Exception {
		if(args.length==0) {
			logger.info("Usage: <program> <input.2bit> [<seq_name> [<start> [<length>]]]");
			logger.info("Resulting fasta data will be written in stdout.");
			return;
		}
		TwoBitParser p = new TwoBitParser(new File(args[0]));
		if(args.length==1) {
			String[] names = p.getSequenceNames();
			for(int i=0;i<names.length;i++) {
				p.setCurrentSequence(names[i]);
				p.printFastaSequence();
				p.close();
			}
		}
		else {
			String name = args[1];
			p.setCurrentSequence(name);
			if(args.length>2) {
				long start = Long.parseLong(args[2]);
				p.skip(start);
			}
			if(args.length>3) {
				long len = Long.parseLong(args[3]);
				p.printFastaSequence(len);
			}
			else {
				p.printFastaSequence();
			}
			p.close();
		}
		p.closeParser();
	}
}

