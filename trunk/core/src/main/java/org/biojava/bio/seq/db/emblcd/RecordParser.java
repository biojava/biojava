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

package org.biojava.bio.seq.db.emblcd;

/**
 * <code>RecordParser</code> is a general purpose parser for arrays of
 * bytes read from any the the EMBL CD-ROM format index files.
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.2
 */
class RecordParser
{
    // File number field (division.lkp & entrynam.idx)
    private byte [] fNumBytes = new byte [2];
    // Offset fields (entrynam.idx)
    private byte [] rPosBytes = new byte [4];
    private byte [] sPosBytes = new byte [4];
    // Record number fields (acnum.trg & acnum.hit)
    private byte [] recNumBytes = new byte [4];
    private byte [] recTotBytes = new byte [4];

    /**
     * Creates a new <code>RecordParser</code> object.
     */
    RecordParser() { }

    /**
     * <code>parseInt4</code> creates a long from Little-endian. Named
     * after the EMBOSS Ajax function which wrote the data.
     *
     * @param int4 a <code>byte []</code> array.
     *
     * @return a <code>long</code>.
     */
    long parseInt4(byte [] int4)
    {
        int result = 0;

        for (int i = 4; --i >= 0;)
        {
            if (int4[i] != 0)
                result += ((int4[i] & 0xff) << (8 * i));
        }

        return result;
    }

    /**
     * <code>parseInt2</code> creates an int from Little-endian. Named
     * after the EMBOSS Ajax function which wrote the data.
     *
     * @param int2 a <code>byte []</code> array.
     *
     * @return an <code>int</code>.
     */
    int parseInt2(byte [] int2)
    {
        int result = 0;

        for (int i = 2; --i >= 0;)
        {
            if (int2[i] != 0)
                result += ((int2[i] & 0xff) << (8 * i));
        }

        return result;
    }

    /**
     * <code>parseDate</code> parses a <code>String</code> from an
     * array of bytes. The date is stored in 4 bytes: 0, unused; 1,
     * year; 2, month; 3, day. With a 1 byte year it's not very much
     * use and I'm not sure that the EMBOSS programs set the value
     * correctly anyway.
     *
     * @param sb a <code>StringBuffer</code>.
     * @param dbDate a <code>byte []</code> array.
     *
     * @return a <code>String</code>.
     */
    String parseDate(StringBuffer sb, byte [] dbDate)
    {
        // The first byte is unused
        for (int i = dbDate.length; --i > 0;)
        {
            sb.append(dbDate[i] + ":");
        }

        // Remove the trailing ':'
        sb.deleteCharAt(sb.length() - 1);
        return sb.substring(0);
    }

    /**
     * <code>parseString</code> parses a <code>String</code> from an
     * array of bytes, skipping the empties.
     *
     * @param sb a <code>StringBuffer</code>.
     * @param characters a <code>byte []</code> array.
     *
     * @return a <code>String</code>.
     */
    String parseString(StringBuffer sb, byte [] characters)
    {
        for (int i = 0; i < characters.length; i++)
        {
            if (characters[i] == 0)
                break;

            sb.append((char) characters[i]);
        }

        return sb.substring(0);
    }

    /**
     * <code>parseDivRecord</code> parses the raw bytes into an array
     * of Objects.
     *
     * @param divRecord a <code>byte []</code> array.
     *
     * @return an <code>Object []</code> array.
     */
    Object [] parseDivRecord(byte [] divRecord)
    {
        // The variable part of record is the name. Other parts are
        // int which sum to 2
        int nameLen = divRecord.length - 2;
        byte [] fNameBytes = new byte [nameLen];

        System.arraycopy(divRecord, 0, fNumBytes, 0, 2);
        System.arraycopy(divRecord, 2, fNameBytes, 0, nameLen);

        StringBuffer sb = new StringBuffer(512);
        Integer fileNumber = new Integer(parseInt2(fNumBytes));
        String    fileName = parseString(sb, fNameBytes);

        return new Object [] { fileNumber, fileName };
    }

    /**
     * <code>parseEntryNamRecord</code> parses the raw bytes into an
     * array of Objects.
     *
     * @param divRecord a <code>byte []</code> array.
     *
     * @return an <code>Object []</code> array.
     */
    Object [] parseEntryNamRecord(byte [] enRecord)
    {
        // The variable part of record is the id. Other parts are
        // long, long, int which sum to 10
        int idLen = enRecord.length - 10;
        byte [] idBytes =  new byte [idLen];

        System.arraycopy(enRecord, 0,           idBytes, 0, idLen);
        System.arraycopy(enRecord, idLen,     rPosBytes, 0, 4);
        System.arraycopy(enRecord, idLen + 4, sPosBytes, 0, 4);
        System.arraycopy(enRecord, idLen + 8, fNumBytes, 0, 2);

        StringBuffer sb = new StringBuffer(512);
        String       seqID = parseString(sb, idBytes);
        Long     rPosition = new Long(parseInt4(rPosBytes));
        Long     sPosition = new Long(parseInt4(sPosBytes));
        Integer fileNumber = new Integer(parseInt2(fNumBytes));

        return new Object [] { seqID, rPosition, sPosition, fileNumber };
    }

    /**
     * <code>parseAcnumTrgRecord</code> parses the raw bytes into an
     * array of Objects.
     *
     * @param divRecord a <code>byte []</code> array.
     *
     * @return an <code>Object []</code> array.
     */
    Object [] parseAcnumTrgRecord(byte [] acRecord)
    {
        // The variable part of record is the acc. Other parts are
        // long, long which sum to 8
        int accLen = acRecord.length - 8;
        byte [] accBytes = new byte [accLen];

        System.arraycopy(acRecord, 0, recNumBytes, 0, 4);
        System.arraycopy(acRecord, 4, recTotBytes, 0, 4);
        System.arraycopy(acRecord, 8, accBytes,    0, accLen);

        StringBuffer sb = new StringBuffer(512);
        Long  rNumber = new Long(parseInt4(recNumBytes));
        Long   rTotal = new Long(parseInt4(recTotBytes));
        String seqAcc = parseString(sb, accBytes);

        return new Object [] { rNumber, rTotal, seqAcc };
    }

    /**
     * <code>parseAcnumHitRecord</code> parses the raw bytes into an
     * array of Objects.
     *
     * @param divRecord a <code>byte []</code> array.
     *
     * @return an <code>Object []</code> array.
     */
    Object [] parseAcnumHitRecord(byte [] rnRecord)
    {
        System.arraycopy(rnRecord, 0, recNumBytes, 0, 4);

        Long rNumber = new Long(parseInt4(recNumBytes));

        return new Object [] { rNumber };
    }
}
