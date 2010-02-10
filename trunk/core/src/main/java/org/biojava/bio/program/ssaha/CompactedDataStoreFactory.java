package org.biojava.bio.program.ssaha;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.nio.IntBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Packing;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.Constants;

/**
 * <p>
 * Builder for a data store that is backed by a java.nio.MappedByteBuffer.
 * This has a limitation that the total size of the mapped buffer and
 * therefore the hash table can not exceed 2 gigs.
 * </p>
 *
 * <p>
 * The data store file has the following structure.
 * <pre>
 * file: header, hash table, nameArray, nameTable, hitTable
 *
 * header:
 *   int hashTablePos, // byte offset in file
 *   int hitTablePos,  // byte offset in file
 *   int nameArrayPos, // byte offset in file
 *   int nameTablePos, // byte offset in file
 *   int wordLength,
 *   int serializedPackingLength,
 *   byte[] serializedPacking
 *
 *   hash table:
 *     int hashTableLength,
 *     int[hashTableLength] hits // index into hitTable
 *
 *  nameArray:
 *    int nameArrayLength,
 *    int[nameArrayLength] nameRecord // byte offset into nameTable
 *
 *  nameRecord:
 *    int nameTableOffset
 *    int sequenceStartOffset
 * 
 *  nameTable:
 *    int nameTableSize, // size in bytes
 *    (short nameLength, char[nameLength] name)[nameTableSize] names
 *
 *  hitTable:
 *    int hitTableSize, // size in bytes
 *    hitTableRecord[hitTableSize] hits
 *
 *  hitTableRecord:
 *    int hitCount,
 *    hitRecord[hitCount] hit
 *
 *  hit:
 *    int offset    // offset into the sequence
 * </pre>
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */

public class CompactedDataStoreFactory implements DataStoreFactory {
  public DataStore getDataStore(File storeFile)
      throws IOException 
  {
      return new CompactedDataStore(storeFile);
  }
  
  public DataStore buildDataStore(
    File storeFile,
    SequenceDB seqDB,
    Packing packing,
    int wordLength,
    int threshold
  ) throws
    IllegalAlphabetException,
    IOException,
    BioException
  {
      return this.buildDataStore(storeFile,
				 new SequenceStreamer.SequenceDBStreamer(seqDB),
				 packing,
				 wordLength,
				 1,
				 threshold);
  }

  public DataStore buildDataStore(
    File storeFile,
    SequenceStreamer streamer,
    Packing packing,
    int wordLength,
    int stepSize,
    int threshold
  ) throws
    IllegalAlphabetException,
    IOException,
    BioException
  { 
    ByteArrayOutputStream packingStream = new ByteArrayOutputStream();
    ObjectOutputStream packingSerializer = new ObjectOutputStream(packingStream);
    packingSerializer.writeObject(packing);
    packingSerializer.flush();
    
    final int structDataSize =
      6 * Constants.BYTES_IN_INT +
      packingStream.toByteArray().length;
    
    final int hashTablePos;
    final int hitTablePos;
    final int nameArrayPos;
    final int nameTablePos;
    
    storeFile.createNewFile();
    final RandomAccessFile store = new RandomAccessFile(storeFile, "rw");
    final FileChannel channel = store.getChannel();
    
    // allocate array for k-tuple -> hit list
    //System.out.println("Word length:\t" + wordLength);
    int words = 1 << (
      (int) packing.wordSize() *
      (int) wordLength
    );
    //System.out.println("Words:\t" + words);
    
    hashTablePos = structDataSize;
    int hashTableSize =
      (int) Constants.BYTES_IN_INT + // hash table length
      words * (int) Constants.BYTES_IN_INT; // hash table entries
    
    //System.out.println("Allocated:\t" + hashTableSize);
    if(hashTableSize < words) {
      throw new AssertionFailure(
        "Possible underflow. number of words: " + words +
        "\tsize of hash table: " + hashTableSize +
        "\tcompared to Integer.MAX_VALUE " + Integer.MAX_VALUE);
    }

    final MappedByteBuffer hashTable_MB = channel.map(
      FileChannel.MapMode.READ_WRITE,
      hashTablePos,
      hashTableSize
    );
    final IntBuffer hashTable = hashTable_MB.asIntBuffer();
    hashTable.put(0, hashTableSize); // write length of k-tuple array
    
    // initialize counts to zero
    for(int i = 0; i < words; i++) {
      hashTable.put(i+1, 0);
    }
    hashTable.position(0);
    
    // System.err.println("And so it begins...");

    // 1st pass
    // writes counts as ints for each k-tuple
    // count up the space required for sequence names
    //

    FirstPassListener fpl = new FirstPassListener(packing, wordLength, stepSize, hashTable);
    streamer.reset();
    while (streamer.hasNext()) {
	streamer.streamNext(fpl);
    }
    
    // map the space for sequence index->name
    //
    nameArrayPos = hashTablePos + hashTableSize;
    int nameArraySize = ((fpl.seqCount * 2) + 1) * Constants.BYTES_IN_INT;
    //System.out.println("seqCount:\t" + seqCount);
    //System.out.println("nameArraySize:\t" + nameArraySize);
    final MappedByteBuffer nameArray_MB = channel.map(
      FileChannel.MapMode.READ_WRITE,
      nameArrayPos,
      nameArraySize
    );
    final IntBuffer nameArray = nameArray_MB.asIntBuffer();
    nameArray.put(0, nameArraySize);
    
    // map the space for sequence names as short length, char* name
    //
    nameTablePos = nameArrayPos + nameArraySize;
    int nameTableSize =
      Constants.BYTES_IN_INT +
      fpl.seqCount * Constants.BYTES_IN_INT +
      fpl.nameChars * Constants.BYTES_IN_CHAR;
    //System.out.println("nameTableSize:\t" + nameTableSize);
    final MappedByteBuffer nameTable = channel.map(
      FileChannel.MapMode.READ_WRITE,
      nameTablePos,
      nameTableSize
    );
    nameTable.putInt(0, nameTableSize);
    nameTable.position(Constants.BYTES_IN_INT);
    
    // add up the number of k-tuples
    //
    int kmersUsed = 0;
    int hitCount = 0;
    for(int i = 0; i < words; i++) {
      int counts = hashTable.get(i + 1);
      if(counts > 0 && counts < threshold) {
        hitCount++;
        kmersUsed += counts;
      }
    }
    
    // map the space for hits
    hitTablePos = nameTablePos + nameTableSize;
    long hitTableSize =
      (long) Constants.BYTES_IN_INT +                            // size
      (long) kmersUsed * (Constants.BYTES_IN_INT) +              // list elements
      (long) hitCount * Constants.BYTES_IN_INT;                  // size of lists
    //System.out.println("hitTableSize:\t" + hitTableSize);
    //System.out.println("hitTableSize:\t" + (int) hitTableSize);
    //System.out.println("hitTablePos:\t" + hitTablePos);
    final MappedByteBuffer hitTable = channel.map(
      FileChannel.MapMode.READ_WRITE,
      hitTablePos,
      (int) hitTableSize
    );
    hitTable.putInt(0, (int) hitTableSize);
    hitTable.position(Constants.BYTES_IN_INT);
    
    // write locations of hit arrays
    int hitOffset = 0;
    for(int i = 0; i < words; i++) {
      int counts = hashTable.get(i+1);
      if(counts > 0 && counts < threshold) {
        try {
        // record location of a block of the form:
        // n,(seqID,offset)1,(seqID,offset)2,...,(seqID,offset)n
        if(hitOffset < 0) {
          throw new IndexOutOfBoundsException("Hit offset negative");
        }
        hashTable.put(i + 1, hitOffset); // wire hash table to hit table
        hitTable.putInt(hitOffset + Constants.BYTES_IN_INT, 0); // initialy we have no hits
        hitOffset +=
          Constants.BYTES_IN_INT +
          counts * (Constants.BYTES_IN_INT);
        } catch (IndexOutOfBoundsException e) {
          System.out.println("counts:\t" + counts);
          System.out.println("word:\t" + i);
          System.out.println("hitOffset:\t" + hitOffset);
          throw e;
        }
      } else if (counts == 0) {
	// nothing - set the number of hits to the flag value -1
        hashTable.put(i + 1, -1);
      } else {
        // too many hits - set the number of hits to the flag value -2
        hashTable.put(i + 1, -2);
      }
    }
    
    // System.err.println("Second pass...");

    // 2nd parse
    // write sequence array and names
    // write hitTable
    
    SecondPassListener spl = new SecondPassListener(packing,
						    wordLength,
						    stepSize,
						    hashTable,
						    nameArray,
						    nameTable,
						    hitTable);
    streamer.reset();
    while (streamer.hasNext()) {
	streamer.streamNext(spl);
    }
    
    //validateNames(seqCount, nameArray, nameTable);
    
    final MappedByteBuffer rootBuffer = channel.map(
      FileChannel.MapMode.READ_WRITE,
      0,
      structDataSize
    );
    
    rootBuffer.position(0);
    rootBuffer.putInt(hashTablePos);
    rootBuffer.putInt(hitTablePos);
    rootBuffer.putInt(nameArrayPos);
    rootBuffer.putInt(nameTablePos);
    rootBuffer.putInt(wordLength);
    rootBuffer.putInt(packingStream.toByteArray().length);
    rootBuffer.put(packingStream.toByteArray());
    
    rootBuffer.force();
    hashTable_MB.force();
    hitTable.force();
    nameArray_MB.force();
    nameTable.force();
    
    return getDataStore(storeFile);
  }
  
    private abstract class PackingListener extends SeqIOAdapter {
	private final Packing packing;
	private final int wordLength;
	private final int stepSize;
	private int pos = -1;
	private int word = 0;
	private int lengthFromUnknown = 0;

	public PackingListener(Packing packing,
			       int wordLength,
			       int stepSize) 
	{
	    this.packing = packing;
	    this.wordLength = wordLength;
	    this.stepSize = stepSize;
	}

	public void startSequence() 
	    throws ParseException
	{
	    pos = 0;
	    word = 0;
	    lengthFromUnknown = 0;
	}

	public void endSequence()
	    throws ParseException
	{
	    foundLength(pos);
	    pos = -1;
	}

	public void foundLength(int length)
	    throws ParseException
	{
	}

	public abstract void processWord(int word, int pos)
	    throws ParseException;

	public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
	    throws IllegalAlphabetException
	{
	    if (alpha != packing.getAlphabet()) {
		throw new IllegalAlphabetException("Alphabet " + alpha.getName() + " doesn't match packing");
	    }

	    int stepCounter = stepSize;
	    for (int i = start; i < (start + length); ++i) {
		word = word >> (int) packing.wordSize();
		try {
		    int p = packing.pack(syms[i]);
		    if (p < 0) {
			lengthFromUnknown = 0;
		    } else {
			lengthFromUnknown++;
			word |= (int) p << ((int) (wordLength - 1) * packing.wordSize());
		    }
		} catch (IllegalSymbolException ex) {
		    throw new BioRuntimeException(ex);
		}

		++pos;
		// System.out.println("Pos = " + pos + "        lengthFromUnknown = " + lengthFromUnknown);
		if (--stepCounter == 0) {
		    stepCounter = stepSize;
		    if (lengthFromUnknown >= wordLength) {
			try {
			    processWord(word, pos - wordLength + 1);
			} catch (ParseException ex) {
			    throw new BioRuntimeException(ex);
			}
		    }
		}
	    }
	}
    }

    private class FirstPassListener extends PackingListener {
	private final IntBuffer hashTable;
        int seqCount = 0;
	int nameChars = 0;

	FirstPassListener(Packing packing,
			  int wordLength,
			  int stepSize,
			  IntBuffer hashTable) 
	{
	    super(packing, wordLength, stepSize);
	    this.hashTable = hashTable;
	}

	public void startSequence()
	    throws ParseException
	{
	    super.startSequence();
	    ++seqCount;
	}

	public void setName(String name) 
	    throws ParseException
	{
            //System.err.println(this + " setting name to " + name);
	    nameChars += name.length();
	}

	public void processWord(int word, int pos)
	    throws ParseException
	{
	    addCount(hashTable, word);
	}
    }

    private class SecondPassListener extends PackingListener {
	private final IntBuffer hashTable;
	private final IntBuffer nameArray;
	private final MappedByteBuffer nameTable;
	private final MappedByteBuffer hitTable;

	private int seqNumber = 0;
	private int concatOffset = 0;

	private String name = ""; // fixme: we need to be cleverer
	private int length = -1;

//  	public void startSequence()
//  	    throws ParseException
//  	{
//  	    super.startSequence();
//  	    System.out.println("Starting sequence");
//  	}

	SecondPassListener(Packing packing, 
			   int wordLength, 
			   int stepSize,
			   IntBuffer hashTable,
			   IntBuffer nameArray,
			   MappedByteBuffer nameTable,
			   MappedByteBuffer hitTable) 
	{
	    super(packing, wordLength, stepSize);

            if( (hashTable == null) ||
                (nameArray == null) ||
                (nameTable == null) ||
                (hitTable  == null) )
            {
              throw new NullPointerException(
                "Buffers must not be null. " +
                "\thashTable: " + hashTable +
                "\tnameArray: " + nameArray +
                "\tnameTable: " + nameTable +
                "\thitTable: " + hitTable );
            }

	    this.hashTable = hashTable;
	    this.nameArray = nameArray;
	    this.nameTable = nameTable;
	    this.hitTable = hitTable;

	    nameTable.position(Constants.BYTES_IN_INT);
	}

	public void setName(String name) {
            //System.err.println(this + " setting name to " + name);
	    this.name = name;
	}

	public void foundLength(int length) {
	    this.length = length;
	}

	public void endSequence() 
	    throws ParseException
	{
	    super.endSequence();
	    
	    nameArray.put((seqNumber * 2) + 1, nameTable.position()-Constants.BYTES_IN_INT);
	    nameArray.put((seqNumber * 2) + 2, concatOffset);
	    // write sequence name length and chars into nameTable
	    nameTable.putInt(name.length());
	    for(int j = 0; j < name.length(); j++) {
		nameTable.putChar((char) name.charAt(j));
	    }

	    ++seqNumber;
	    concatOffset += (length + 100);
	}

	public void processWord(int word, int pos)
	    throws ParseException
	{
	    if (pos < 1) {
		throw new ParseException("pos < 1");
	    }
	    writeRecord(hashTable, hitTable, pos + concatOffset, seqNumber, word);
	}
    }

  private void addCount(IntBuffer buffer, int word) {
    int count = buffer.get(word+1);
    count++;
    buffer.put(word+1, count);
  }
  
  private void writeRecord(
    IntBuffer hashTable,
    MappedByteBuffer hitTable,
    int offset,
    int seqNumber,
    int word
  ) {
    int kmerPointer = hashTable.get(word+1);
    if(kmerPointer >= 0) {
      kmerPointer += Constants.BYTES_IN_INT;

      int hitCount = hitTable.getInt(kmerPointer);
      int pos = kmerPointer + hitCount * (Constants.BYTES_IN_INT) + Constants.BYTES_IN_INT;
      
      hitTable.position(pos);
      hitTable.putInt(offset);
      hitTable.putInt(kmerPointer, hitCount + 1);
    }
  }
  
}
