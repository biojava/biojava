package org.biojava.utils.io;

import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import org.biojava.utils.Constants;

/**
 * Wrapper arround MappedByteBuffers to allow long-indexed access to files
 * larger than 2 gigs.
 *
 * @author Matthews Pocock
 */
public class LargeBuffer {
  /*
   * We will set up MappedByteBuffers that are responsible for PAGE_SIZE
   * bytes. Unfortunately, word boundaries are not aligned, so someone could
   * try to write a double to the last byte in a buffer. So, 
   */

  private static final long PAGE_SIZE;
  private static final long PAGE_OVERLAP;
  
  static {
    PAGE_OVERLAP = Constants.BYTES_IN_LONG;
    PAGE_SIZE = Integer.MAX_VALUE / 8 - PAGE_OVERLAP;
  }
  
  private final long pos;
  private final long size;
  private final FileChannel channel;
  private final FileChannel.MapMode mode;

  private long position = 0;
  private int lastBufferIndex = -1;
  private MappedByteBuffer lastBuffer = null;
  
  
  public LargeBuffer(
    FileChannel channel,
    FileChannel.MapMode mode,
    long pos,
    long size
  ) throws IOException {
    this.channel = channel;
    this.mode = mode;
    this.pos = pos;
    this.size = size;
  }
  
  private MappedByteBuffer getBuffer(int index)
  throws IOException {
    if(index != lastBufferIndex) {
      System.out.println("Allocating page: " + index);
      long offset = PAGE_SIZE * index;
      System.out.println("From: " + (pos + offset));
      System.out.println("Size: " + Math.min(size - offset, PAGE_SIZE + PAGE_OVERLAP));
      lastBuffer = channel.map(
        mode,
        pos + offset,
        Math.min(size - offset, PAGE_SIZE + PAGE_OVERLAP)
      );
      lastBufferIndex = index;
      System.out.println("Done");
    }
    
    return lastBuffer;
  }
  
  public byte get(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.get(offset);
  }
  
  public byte get()
  throws IndexOutOfBoundsException, IOException {
    byte val = get(position);
    position += Constants.BYTES_IN_BYTE;
    return val;
  }
  
  public void put(long pos, byte b)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.put(offset, b);
  }
  
  public void put(byte val)
  throws IndexOutOfBoundsException, IOException {
    put(position, val);
    position += Constants.BYTES_IN_BYTE;
  }
  
  public char getChar(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getChar(offset);
  }
  
  public char getChar()
  throws IndexOutOfBoundsException, IOException {
    char val = getChar(position);
    position += Constants.BYTES_IN_CHAR;
    return val;
  }
  
  public void putChar(long pos, char c)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putChar(offset, c);
  }
  
  public void putChar(char val)
  throws IndexOutOfBoundsException, IOException {
    putChar(position, val);
    position += Constants.BYTES_IN_CHAR;
  }
  
  public double getDouble(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getDouble(offset);
  }
  
  public double getDouble()
  throws IndexOutOfBoundsException, IOException {
    double val = getDouble(position);
    position += Constants.BYTES_IN_DOUBLE;
    return val;
  }
  
  public void putDouble(long pos, double d)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putDouble(offset, d);
  }
  
  public void putDouble(double val)
  throws IndexOutOfBoundsException, IOException {
    putDouble(position, val);
    position += Constants.BYTES_IN_DOUBLE;
  }
  
  public float getFloat(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getFloat(offset);
  }
  
  public float getFloat()
  throws IndexOutOfBoundsException, IOException {
    float val = getFloat(position);
    position += Constants.BYTES_IN_FLOAT;
    return val;
  }
  
  public void putFloat(long pos, float f)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putFloat(offset, f);
  }
  
  public void putFloat(float val)
  throws IndexOutOfBoundsException, IOException {
    putFloat(position, val);
    position += Constants.BYTES_IN_FLOAT;
  }
  
  public int getInt(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getInt(offset);
  }
  
  public int getInt()
  throws IndexOutOfBoundsException, IOException {
    int val = getInt(position);
    position += Constants.BYTES_IN_INT;
    return val;
  }
  
  public void putInt(long pos, int i)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putInt(offset, i);
  }
  
  public void putInt(int val)
  throws IndexOutOfBoundsException, IOException {
    putInt(position, val);
    position += Constants.BYTES_IN_INT;
  }
  
  public long getLong(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getLong(offset);
  }
  
  public long getLong()
  throws IndexOutOfBoundsException, IOException {
    long val = getLong(position);
    position += Constants.BYTES_IN_LONG;
    return val;
  }
  
  public void putLong(long pos, long l)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putLong(offset, l);
  }
  
  public void putLong(long val)
  throws IndexOutOfBoundsException, IOException {
    putLong(position, val);
    position += Constants.BYTES_IN_LONG;
  }
  
  public short getShort(long pos)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    return buffer.getShort(offset);
  }
  
  public short getShort()
  throws IndexOutOfBoundsException, IOException {
    short val = getShort(position);
    position += Constants.BYTES_IN_SHORT;
    return val;
  }
  
  public void putShort(long pos, short s)
  throws IndexOutOfBoundsException, IOException {
    int offset = getOffset(pos);
    int index = getIndex(pos);
    
    MappedByteBuffer buffer = getBuffer(index);
    buffer.putShort(offset, s);
  }
  
  public void putShort(short val)
  throws IndexOutOfBoundsException, IOException {
    putShort(position, val);
    position += Constants.BYTES_IN_SHORT;
  }
  
  public long position() {
    return position;
  }
  
  public void position(long pos) {
    this.position = pos;
  }
  
  private int getOffset(long pos)
  throws IndexOutOfBoundsException {
    if(pos > size) {
      throw new IndexOutOfBoundsException();
    }
    return (int) (pos % PAGE_SIZE);
  }
  
  private int getIndex(long pos) {
    return (int) (pos / (long) PAGE_SIZE);
  }
  
  public void force() {
//    for(Iterator i = buffers.iterator(); i.hasNext(); ) {
//      MappedByteBuffer buff = (MappedByteBuffer) i.next();
//      buff.force();
//    }
  }
}
