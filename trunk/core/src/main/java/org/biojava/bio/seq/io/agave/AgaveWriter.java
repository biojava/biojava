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
package org.biojava.bio.seq.io.agave;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Iterator;

import org.biojava.bio.Annotatable;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleAssembly;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;

/**
 * Writes Sequence into AGAVE XML document.  The AGAVE
 * format is defined in agave.dtd which can be downloaded
 * from http://www.agavexml.org.
 *
 * @author Hanning Ni
 * @author Brian King
 */
public class AgaveWriter
{
    /**
     * Implements indenting for elements.
     * @author Brian King
     */
    protected class Indent
    {
        /** Base indent    */
        protected String mIndent = new String(INDENT);

        /**
         * Add a level of indentation
        */
        public void
        indent()
        {
            mIndent += INDENT;
        }

        /**
         * Remove a level of indentation
         *
         */
        public void
        unIndent()
        {
            int len  = mIndent.length();
            int iLen = INDENT.length();
            if(len >= iLen)
            {
                mIndent = mIndent.substring(0, len - iLen);
            }
        }

        /**
         * Return the current indent
         *
         */
        public String
        toString()
        {
            return mIndent;
        }
    }

   /** use a two space indent */
    public static final String INDENT = "  ";

    /** Write to XML document    */
    protected PrintWriter mOut;

    /** indent */
    protected Indent mIndent;

    /** writes PCDATA replacing XML characters with escape entities */
    protected PCDATAFilterWriter mFilter;

    /**
     * The AnnotationMap to use for getting
     * AGAVE XML attributes from a Sequence Annotation.
     */
    protected AGAVEAnnotFilter mAnnotFilter;

    /** write DOCTYPE if true */
    protected boolean mWriteDocType = true;

    /**
     * Default constructor uses generic annotation to attribute mapping.
     *
     */
    public AgaveWriter()
    {
        mAnnotFilter = SimpleAnnotFilter.SIMPLE_ANNOT_FILTER_FACTORY.getInstance();
        mIndent = new Indent();
    }

    /**
     * Construct with data source specific annotation to attribute
     * mapping.
     *
     */
    public AgaveWriter(AGAVEAnnotFilter filter)
    {
        mAnnotFilter = filter;
        mIndent = new Indent();
    }

    /**
     * Set flag that determines if XML DOCTYPE is written
     * or not.  Default is true.
     *
     */
    public void
    setWriteDocType(boolean writeDocType)
    {
        mWriteDocType = writeDocType;
    }

     /**
     * Write sequence into AGAVE XML format.
     * @param seq maybe the  or simple sequence
     * <pre>
     * if annotation of seq has chromosome information , generate <chromosome> tag
     * if seq is SimpleAssembly, generate <contig> tag
     * otherwise, generate <bio_sequence> tag
     * currently  each top-level sequence is wrapped into seperated sciobj xml file
     * </pre>
     */
    public void writeSequence(Sequence seq, PrintStream os)
           throws IOException
    {
        mOut = new PrintWriter(os);
        mFilter = new PCDATAFilterWriter(mOut);

        mOut.println("<?xml version=\"1.0\"?>");
        if (mWriteDocType)
        {
            mOut.println("<!DOCTYPE sciobj SYSTEM \"agave.dtd\">");
        }
        writeHeader();
        write(seq);
        writeFooter();
    }

    /**
     * Write &lt;sciobj&gt;
     *
     */
    protected void
    writeHeader()
    {
        mOut.println("<sciobj version=\"2\">");
    }

    /**
     * Write &lt;/sciobj&gt;
     */
    protected void
    writeFooter()
    {
        mOut.println("</sciobj>");
        mOut.flush();
    }

    /**
     *
     * Writing Sequence.
     * @param seq is simple sequence or simple assembly
     *
     */
    protected void
    write(Sequence seq) throws IOException
    {
        String chrom_num = mAnnotFilter.getChromNum( seq.getAnnotation() );
        if( chrom_num != null )
        {
            mIndent.indent();
            mOut.print(mIndent);
            mOut.println("<chromosome number=\"" + chrom_num + "\">");
            writeContig((Annotatable) seq );
            mOut.print(mIndent);
            mOut.println("</chromosome>");
            mIndent.unIndent();
        }
        else
        {
            if( seq instanceof SimpleAssembly)
            {
               writeContig(  (Annotatable)seq ) ;
            }
            else
            {
                writeBioSequence((Annotatable)seq);
            }
        }
    }

    /**
     *
     *
     */
    protected void
    writeContig(Annotatable seq ) throws IOException
    {
    /**
    <!ELEMENT contig (db_id , view? , note? , fragment_order* , unordered_fragments? ,  assembly? ,
                  sequence? , sequence_map* ,  map_location* )>
    <!ATTLIST contig length NMTOKEN  #REQUIRED >
    */
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<contig");


        mOut.print(" length=\"");
        mOut.print( ((Sequence)seq).length() );
        mOut.print('"');
        mOut.println(">");

        writeDbId(seq);
        writeNote(seq);
        writeAssembly(seq);
        writeDNA(seq) ;
        writeMapLocation(seq) ;


        mOut.print(mIndent);
        mOut.println("</contig>");
        mIndent.unIndent();

    }

    /**
     *
     *
     */
    protected void
    writeAssembly(Annotatable seq) throws IOException
    {
        /** <!ELEMENT assembly (bio_sequence | fragment_order)+ > */
        mIndent.indent();
        mOut.print(mIndent);
        mOut.println("<assembly>");
        if(seq instanceof SimpleSequence)
            write((Sequence) seq ) ;
        else if( seq instanceof SimpleAssembly)
        {
            for (Iterator i  =((Sequence) seq).features(); i.hasNext(); )
            {
                Feature subf = (Feature)i.next();
                if( subf instanceof ComponentFeature)
                    writeBioSequence((Annotatable)  ((ComponentFeature)subf).getComponentSequence() ) ;
            }
        }

        mOut.print(mIndent);
        mOut.println("</assembly>");
        mIndent.unIndent();
    }

    /**
     *
     *
     */
    protected void
    writeBioSequence(Annotatable seq ) throws IOException
    {
        String s;
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<bio_sequence");


        // write attributes
        //
        if( seq instanceof Sequence)
        {
            mOut.print(" seq_length=\"");
            mOut.print( ((Sequence)seq).length() );
            mOut.print('"');
        }

        s = mAnnotFilter.getOrganism(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" organism_name=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getMolType( seq.getAnnotation() );
	if (s == null && (seq instanceof Sequence)) {
	    s = ((Sequence) seq).getAlphabet().getName();
	}
        if (s != null)
        {
            mOut.print(" molecule_type=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getElementId(seq.getAnnotation());
	if (s == null && (seq instanceof Sequence)) {
	    s = ((Sequence) seq).getName();
	}
        if (s != null)
        {
            mOut.print(" element_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getSequenceId(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" sequence_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getTaxonId(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" taxon_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getCloneId(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" clone_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getCloneLibrary(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" clone_library=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getChromosome(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" chromosome=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getMapPosition(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" map_position=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getEcNumber(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" ec_number=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getCreateDate(seq.getAnnotation());
        if (s != null)
        {
            mOut.print(" update_date=\"");
            mFilter.write(s);
            mOut.print('"');
        }


        mOut.println(">");  // close bio_sequence tag

        // write content
        writeDbId(seq) ;
        writeDNA( seq ) ;
        writeAltIds( seq ) ;
        writeXrefs( seq ) ;
        writeSequenceMap2(seq);   // THOMASD hacked this....
        writeMapLocation(seq) ;
        writeClassification(seq);

        mIndent.unIndent();
        mOut.print(mIndent);
        mOut.println("</bio_sequence>");
    }

    /**
     * group sequence_map by getSource()
     */
    protected void
    writeSequenceMap(Annotatable seq) throws IOException
    {
        for (Iterator i  = ((FeatureHolder)seq).features(); i.hasNext();)
        {
            Feature f = (Feature) i.next();
            String type = f.getSource();
            if( type.equalsIgnoreCase("classification") )
                continue ;
            // writeFeature( f );
            writeSequenceMap2( f );
         }
    }

    /**
     *
     *
     */
    protected void
    writeClassification(Annotatable seq) throws IOException
    {
        for (Iterator i  = ((FeatureHolder)seq).features(); i.hasNext();)
        {
            Feature f = (Feature) i.next();
            String type = f.getSource();
            if( type.equalsIgnoreCase("classification") )
            {
                writeClassification2(f) ;
            }
        }
     }

    /**
     *
     *
     */
     private void
     writeClassification2(Annotatable f) throws IOException
     {
     /**
     <!ELEMENT classification (description? , id_alias* , evidence? )>
     <!ATTLIST classification system      CDATA  #REQUIRED
                         id        CDATA  #REQUIRED
                         type        CDATA  #REQUIRED
                         assigned_by CDATA  #IMPLIED >
     */
        mOut.print(mIndent);
        mOut.print("<classification");
        mIndent.indent();

        // write attributes
        //
        String s = mAnnotFilter.getClassifySystem( f.getAnnotation() ); ;
        if(s != null)
        {
            mOut.print(" system=\"");
            mFilter.write(s);
            mOut.print('"');
        }
        s = mAnnotFilter.getClassifyId( f.getAnnotation() ); ;
        if(s != null)
        {
            mOut.print(" id=\"");
            mFilter.write(s);
            mOut.print('"');
        }
        s = mAnnotFilter.getClassifyType( f.getAnnotation() ); ;
        if(s != null)
        {
            mOut.print(" type=\"");
            mFilter.write(s);
            mOut.print('"');
        }
        mOut.println(">");  // close map tag

        writeDescription( f ) ;
        writeIdAlias( f ) ;
        writeEvidence( f ) ;

        mOut.print(mIndent);
        mOut.println("</classification>");
        mIndent.unIndent();
     }

     /**
      *
      *
      */
     private void
     writeIdAlias(Annotatable f) throws IOException
     {
       AGAVEIdAlias[] annots =  mAnnotFilter.getIdAlias(f.getAnnotation());
       if( annots != null )
       {
           mIndent.indent();
           for(int i = 0 ; i < annots.length; i++)
              mOut.print( annots[i].toString(mIndent.toString(), INDENT ));
           mIndent.unIndent();
       }
     }

     /**
     *  Write SequenceMap XML
     *
     */
    protected void
    writeSequenceMap2(Annotatable f) throws IOException
    {
/*
<!ELEMENT sequence_map  (note? , computation? , annotations? )>
<!ATTLIST sequence_map label        CDATA #IMPLIED >
*/
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<sequence_map");


        // write attributes
        //
        String label = mAnnotFilter.getLabel(f.getAnnotation()) ;
        if(label != null)
        {
            mOut.print(" label=\"");
            mFilter.write(label);
            mOut.print('"');
        }
        mOut.println(">");  // close map tag

        // write content
        //

        //write note
        //
        writeNote(f) ;

        //ignore write computation

        writeAnnotations((FeatureHolder) f);

        mOut.print(mIndent);
        mOut.println("</sequence_map>");
        mIndent.unIndent();

    }

    /**
     *
     *
     */
    protected void writeAnnotations(FeatureHolder f) throws IOException   // changed to FeatureHolder THOMASD
    {
     //write annotation
        //    <!ELEMENT annotations       (seq_feature | gene | comp_result )+ >
        mIndent.indent();
        mOut.print(mIndent);
        mOut.println("<annotations>");

  /*      for (Iterator i  = f.features(); i.hasNext(); )
        {
                Feature feature = (Feature) i.next();
                if( feature.getSource().equalsIgnoreCase("annotations") )
                {
                    f = feature  ;
                    break;
                }
        } */
        for (Iterator i  = f.features(); i.hasNext(); )
        {
                Feature feature = (Feature) i.next();
                String type = feature.getSource() ;
                if( type == null || type.equalsIgnoreCase( "seq_feature" ) )
                   writeSeqFeature(feature);
                else if( type.equalsIgnoreCase( "gene") )
                   writeGene( feature);
                else { //default mapping to comp_result
		    // writeCompResult( feature);
		   writeSeqFeature(feature);  // THOMASD
		}
        }
        mOut.print(mIndent);
        mOut.println("</annotations>");
        mIndent.unIndent();
    }

    /**
     *
     *
     */
   protected void
   writeGene(Annotatable f) throws IOException
   {
   /**
   <!ELEMENT gene (classification* , note? , seq_location , xrefs? ,
                evidence? , qualifier* , seq_feature* , related_annot* ,
                transcript*)>
   <!ATTLIST gene  element_id   ID     #IMPLIED
                label        CDATA  #IMPLIED >
   **/
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<gene");

        String s ;

        s = mAnnotFilter.getElementId( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" element_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getLabel( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" label=\"");
            mFilter.write(s);
            mOut.print('"');
        }


        mOut.println('>');

        // write content
        writeNote(f) ;
        writeSeqLocation( f ) ;
        writeXrefs( f ) ;
        writeEvidence(f) ;
        writeProperty(f,  AGAVEProperty.QUALIFIER ) ;
        writeSubSeqFeature( f ) ;
        writeRelatedAnnot( f ) ;
        writeTranscript( f ) ;


       mOut.print(mIndent);
       mOut.println("</gene>");
       mIndent.unIndent();
   }

   /**
    *
    *
    */
   protected void
   writeTranscript( Annotatable f) throws IOException
   {
      /**
       * <!ELEMENT transcript (exons , cds? , mrna? , predicted_protein?)>
       */
      for (Iterator i  = ((FeatureHolder)f).features(); i.hasNext(); )
      {
          Feature subf = (Feature)i.next();
          if( subf.getSource().equalsIgnoreCase("transcript") )
              writeTranscript2(subf);
      }
   }

   /**
    *
    *
    */
   private void
   writeTranscript2(Annotatable f) throws IOException
   {
         mIndent.indent();
         mOut.print(mIndent);
         mOut.println("<transcript>");
         writeExons(f);
         for (Iterator i  = ((FeatureHolder)f).features(); i.hasNext(); )
         {
            Feature subf = (Feature)i.next();
            String type = subf.getSource() ;
            if( type.equalsIgnoreCase("cds") )
               writeCds(subf);
            else if( type.equalsIgnoreCase("mrna") )
               writeMrna(subf);
            else if( type.equalsIgnoreCase("predicted_protein") )
               writePredictedProtein(subf);

         }
         mOut.print(mIndent);
         mOut.println("</transcript>");
         mIndent.unIndent();
   }

   /**
    *
    *
    */
   private void
   writeExons(Annotatable f) throws IOException
   {
       String[] ids =  mAnnotFilter.getExonIds(f.getAnnotation());
       if( ids != null )
       {
           mIndent.indent();
           mOut.print(mIndent);
           mOut.println("<exons>");
           for(int i = 0 ; i < ids.length; i++)
           {
               mOut.println(mIndent + INDENT + "<element_id id=\"" + ids[i] + "\"/>");
           }
           mOut.print(mIndent);
           mOut.println("</exons>") ;
           mIndent.unIndent();
       }
   }

   /**
    *
    *
    */
   private void writeCds(Annotatable f) throws IOException
   {
       mIndent.indent();
       mOut.print(mIndent);
       mOut.println("<cds>");
       writeBioSequence( f) ;
       mOut.print(mIndent);
       mOut.println("</cds>");
       mIndent.unIndent();
   }

   /**
    *
    *
    */
   private void writeMrna(Annotatable f) throws IOException
   {
       mIndent.indent();
       mOut.print(mIndent);
       mOut.println("<mrna>");
       writeBioSequence( f) ;
       mOut.print(mIndent);
       mOut.println("</mrna>");
       mIndent.unIndent();
   }

   /**
    *
    *
    */
   private void writePredictedProtein(Annotatable f) throws IOException
   {
       mIndent.indent();
       mOut.print(mIndent);
       mOut.println("<predicted_protein>");
       writeBioSequence( f) ;
       mOut.print(mIndent);
       mOut.println("</predicted_protein>");
       mIndent.unIndent();
   }

    /**
     *  Write SeqFeature XML
     */
    protected void
    writeSeqFeature( Annotatable f) throws IOException
    {
/*
<!ELEMENT seq_feature  (classification* , note? , seq_location , xrefs? ,
                        evidence? , qualifier* ,  seq_feature* , related_annot*)>
<!ATTLIST seq_feature element_id   ID     #IMPLIED
                      feature_type CDATA  #REQUIRED
                      label        CDATA  #IMPLIED >
*/
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<seq_feature");


        // write attributes
        //
        String s = mAnnotFilter.getFeatureType( f.getAnnotation() );
	if (s == null) {  // THOMASD
	    s = ((Feature) f).getType();
	}
        mOut.print(" feature_type=\"");
        mFilter.write( s == null? "default" : s );
        mOut.print('"');

        s = mAnnotFilter.getElementId( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" element_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getLabel( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" label=\"");
            mFilter.write(s);
            mOut.print('"');
        }


        mOut.println('>');

        // write content
        writeNote(f) ;
        writeSeqLocation( f ) ;
        writeXrefs( f ) ;
        writeEvidence(f) ;
        writeProperty(f,  AGAVEProperty.QUALIFIER ) ;
        writeSubSeqFeature( f ) ;
        writeRelatedAnnot( f ) ;


       mOut.print(mIndent);
       mOut.println("</seq_feature>");
       mIndent.unIndent();

    }

    /**
     *
     *
     */
    private void
    writeRelatedAnnot( Annotatable f)
    {
       AGAVERelatedAnnot[] annots =  mAnnotFilter.getRelatedAnnot(f.getAnnotation());
       if( annots != null )
       {
           for(int i = 0 ; i < annots.length; i++)
              mOut.print( annots[i].toString(mIndent + INDENT, INDENT));
       }
    }

    /**
     *
     *
     */
    private void
    writeEvidence(Annotatable f) throws IOException
    {
        for (Iterator i  = ((FeatureHolder)f).features(); i.hasNext(); )
        {
            Feature subf = (Feature)i.next();
            if( subf.getSource().equalsIgnoreCase("evidence") )
                writeEvidence2(subf);
        }
    }

    /**
     *
     *
     */
    private void
    writeEvidence2(Annotatable f) throws IOException
    {
        //element_ids
        writeElementIds( f ) ;
        //comp_result
        for (Iterator i  = ((FeatureHolder)f).features(); i.hasNext(); )
        {
            Feature subf = (Feature)i.next();
            writeCompResult(  subf );
        }
    }

    /**
     *
     *
     */
    private void
    writeElementIds(Annotatable f) throws IOException
    {
       String[] ids =  mAnnotFilter.getElementIds(f.getAnnotation());
       if( ids != null )
       {
           for(int i = 0 ; i < ids.length; i++)
           {
                mOut.println(mIndent + INDENT + "<element_id id=\"" + ids[i] + "\"/>");
           }
       }
    }

    /*
    <!ELEMENT comp_result  (note? , match_desc? , match_align? , query_region? ,
                        match_region? , result_property* , result_group* ,
                        related_annot*)>
    <!ATTLIST comp_result  element_id           ID       #IMPLIED
                       result_id            NMTOKEN  #IMPLIED
                       group_order          NMTOKEN  #IMPLIED
                       result_type          CDATA    #REQUIRED
                       feature_type         CDATA    #IMPLIED
                       on_complement_strand  (true | false )  'false'
                       confidence           NMTOKEN  #IMPLIED
                       align_length         NMTOKEN  #IMPLIED
                       align_units          (bp | AA) #IMPLIED >
    */
    protected
    void writeCompResult(Annotatable f) throws IOException
    {
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<comp_result");


        // write attributes
        //
        String s = mAnnotFilter.getResultType( f.getAnnotation() );
        mOut.print(" result_type=\"");
        mFilter.write( s == null? "default" : s );
        mOut.print('"');

        s = mAnnotFilter.getElementId( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" element_id=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        String strand = "false";
        if( f instanceof StrandedFeature )
        {
            char mark = ((StrandedFeature)f).getStrand().getToken() ;
            if ( mark == '-' )
                strand= "true" ;
        }
        mOut.print(" on_complement_strand=\"");
        mFilter.write(strand);
        mOut.print('"');

        s = mAnnotFilter.getGroupOrder( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" group_order=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getFeatureType( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" feature_type=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getConfidence( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" confidence=\"");
            mFilter.write(s);
            mOut.print('"');
        }

        s = mAnnotFilter.getAlignLength( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" align_length=\"");
            mFilter.write(s);
            mOut.print('"');
        }
        s = mAnnotFilter.getAlignUnits( f.getAnnotation() );
        if (s != null)
        {
            mOut.print(" align_units=\"");
            mFilter.write(s);
            mOut.print('"');
        }
        mOut.println(">");

       writeNote( f ) ;
       writeMatchDesc( f ) ;
       writeMatchAlign( f ) ;
       writeQueryRegion( f ) ;
       writeMatchRegion( f ) ;
       writeProperty( f, AGAVEProperty.RESULT_PROPERTY) ;
       writeResultGroup( f )  ;
       writeRelatedAnnot(f);

       mOut.print(mIndent);
       mOut.println("</comp_result>");
       mIndent.unIndent();
    }

    /**
     *   <!ELEMENT result_group  (comp_result+ )>
     *   <!ATTLIST result_group  group_order NMTOKEN  '0' >
     */
    private void
    writeResultGroup(Annotatable f) throws IOException
    {
        Iterator i = ((FeatureHolder)f).features() ;
        if( i == null ) //void
            return ;
        while( i.hasNext() )
        {
            Feature subf = (Feature)i.next();
            String type = subf.getSource();
            if( type != null && type.equalsIgnoreCase("result_group") )
            {
                mIndent.indent();
                mOut.println( mIndent + "<result_group>");
                for(Iterator k = subf.features(); k.hasNext();)
                    writeCompResult( (Feature)k.next() ) ;
                mOut.println(mIndent + "</result_group>");
                mIndent.unIndent();
            }
            else
            {
                 mIndent.indent();
                 mOut.println(mIndent + "<result_group>");
                 writeCompResult( subf ) ;
                 mOut.println(mIndent + "</result_group>");
                 mIndent.unIndent();
            }
        }
    }

    /**
     *
     *
     */
    private void
    writeQueryRegion(Annotatable f) throws IOException
    {
        AGAVEQueryRegion s =   mAnnotFilter.getQueryRegion( f.getAnnotation() );
        if( s != null )
        {
            mIndent.indent();
            mOut.println( s.toString(mIndent.toString()  , INDENT ) );
            mIndent.unIndent();
        }
    }

    /**
     *
     *
     */
    private void
    writeMatchRegion(Annotatable f) throws IOException
    {
        AGAVEMatchRegion s =   mAnnotFilter.getMatchRegion( f.getAnnotation() );
        if( s != null )
        {
            mIndent.indent();
            mOut.println( s.toString(mIndent.toString()  , INDENT ) );
            mIndent.unIndent();
        }
    }

    /**
     *
     *
     */
    private void
    writeMatchAlign(Annotatable f) throws IOException
    {
        String s =   mAnnotFilter.getMatchAlign( f.getAnnotation() );
        if( s != null )
        {
            if (s.startsWith("<match_align") )
            {
                mFilter.write(s);
             }else{
                mIndent.indent();
                mOut.print(mIndent);
                mOut.print("<match_align>");
                mFilter.write(s);
                mOut.println("</match_align>");
                mIndent.unIndent();
             }
         }
    }

    /**
     *
     *
     */
   private void
   writeMatchDesc(Annotatable f) throws IOException
    {
        String s =   mAnnotFilter.getMatchDesc( f.getAnnotation() );
        if( s != null )
        {
            if (s.startsWith("<match_desc") )
            {
                mFilter.write(s);
             }else{
                mIndent.indent();
                mOut.print(mIndent);
                mOut.print("<match_desc>");
                mFilter.write(s);
                mOut.println("</match_desc>");
                mIndent.unIndent();
             }
         }
    }

   /**
    *
    *
    */
    private void
    writeSubSeqFeature(Annotatable f ) throws IOException
    {
        for (Iterator i  = ((FeatureHolder)f).features(); i.hasNext(); )
        {
            Feature subf = (Feature)i.next();
            //if the feature is not declared as <evidence>
            //treated as seq_feature
            if( subf.getSource().equalsIgnoreCase("evidence") )
                continue ;
           if( subf.getSource().equalsIgnoreCase("transcript") )
                continue ;
            writeSeqFeature( subf);
        }
    }

    /**
     *
     *
     */
    private void
    writeSeqLocation( Annotatable f) throws IOException
    {
        int min = ((Feature)f).getLocation().getMin();
        int max = ((Feature)f).getLocation().getMax();
        boolean on_complement_strand  = false;
        if( f instanceof StrandedFeature)
        {
            int type =  ((StrandedFeature)f).getStrand().getValue() ;
            if( type == 1 )
                on_complement_strand = false ;
            else if( type == -1 )
                on_complement_strand = true ;
        }
        mIndent.indent();
        mOut.print(mIndent);
        mOut.print("<seq_location");
        mOut.print(" least_start=\"");
        mOut.print(min);
        mOut.print('"');

        mOut.print(" greatest_end=\"");
        mOut.print(max);
        mOut.print('"');

        if( f instanceof StrandedFeature)
        {
           mOut.print(" on_complement_strand=\"");
           mOut.print(on_complement_strand ? "true" : "false");
           mOut.print('"');
        }
        mOut.print(">");

        mOut.print(min + ".." + max);
        mOut.println("</seq_location>");
        mIndent.unIndent();
    }

    /**
     *
     *
     */
    private void
    writeDNA(Annotatable seq) throws IOException
    {
       if( seq instanceof Sequence)
       {
           int i = ((Sequence)seq).length();
           if (i > 0)
          {
              mIndent.indent();
              mOut.print(mIndent);
              mOut.print("<sequence>");
              int j ;
              for (j = 0 ; j < i / 80; j ++)
              {
               mOut.print(((Sequence)seq).subList(80 * j + 1, 80 * (j+1)).seqString() + "\n");
              }

              mOut.print(((Sequence)seq).subList(80 * j + 1, i ).seqString());
              mOut.println("</sequence>");
              mIndent.unIndent();
           }
        }
    }

    /**
     *
     *
     */
    private void writeXrefs(Annotatable f)
    {
       AGAVEXrefs[] xrefs =  mAnnotFilter.getXrefs(f.getAnnotation());
       if( xrefs != null )
       {
           mIndent.indent();
           for(int i = 0 ; i < xrefs.length; i++)
              mOut.print( xrefs[i].toString(mIndent.toString(), INDENT));
           mIndent.unIndent();
       }
    }

    /**
     *
     *
     */
    private void
    writeAltIds(Annotatable f) throws IOException
    {
       AGAVEDbId[] db_id =  mAnnotFilter.getAltIds(f.getAnnotation());
       if( db_id!= null )
       {
           mIndent.indent();
           mOut.println(mIndent + "<alt_ids>");
           mIndent.indent();
           for(int i = 0 ; i < db_id.length; i++)
              mOut.print( db_id[i].toString(mIndent.toString(), INDENT));
           mIndent.unIndent();
           mOut.println(mIndent + "</alt_ids>");
           mIndent.unIndent();
       }
    }

    /**
     *
     *
     */
    private void
    writeProperty(Annotatable f, String type)
    {
        AGAVEProperty[] aps =  mAnnotFilter.getProperty( f.getAnnotation(),type) ;
        if( aps != null )
        {
            mIndent.indent();
            for( int index  = 0 ; index < aps.length; index++)
                mOut.print( aps[index].toString(mIndent.toString(), INDENT));
            mIndent.unIndent();
        }
    }

    /**
     *
     *
     */
    private void
    writeMapLocation(Annotatable f) throws IOException
    {
       AGAVEMapLocation[] mls =  mAnnotFilter.getMapLocation(f.getAnnotation());
       if( mls != null  )
       {
           mIndent.indent();
           for(int i = 0 ; i < mls.length; i++)
              mOut.print( mls[i].toString(mIndent.toString(),INDENT));
           mIndent.unIndent();
       }
    }

    /**
     *
     *
     */
    private void
    writeDbId(Annotatable f) throws IOException
    {
        AGAVEDbId db_id = (AGAVEDbId)mAnnotFilter.getDbId(f.getAnnotation());
        if( db_id!= null )
        {
            mIndent.indent();
            mOut.println( db_id.toString(mIndent.toString(), INDENT) );
            mIndent.unIndent();
        }
    }

    /**
     *
     *
     */
    private void
    writeDescription( Annotatable f) throws IOException
    {
        String s =   mAnnotFilter.getDescription( f.getAnnotation() );
        if( s != null )
        {
            if (s.startsWith("<description") )
            {
                mFilter.write(s);
             }else{
                mIndent.indent();
                mOut.print(mIndent);
                mOut.print("<description>");
                mFilter.write(s);
                mOut.println("</description>");
                mIndent.unIndent();
             }
         }
    }

    /**
     *
     *
     */
    private void
    writeNote( Annotatable f) throws IOException
    {
        String note = mAnnotFilter.getNote(f.getAnnotation());
        if( note!= null )
        {
            if (note.startsWith("<note"))
            {
                mFilter.write(note);
            }
            else
            {
                mIndent.indent();
                mOut.print(mIndent);
                mOut.print("<note>");
                mFilter.write(note);
                mOut.println("</note>");
                mIndent.unIndent();
            }
        }
    }
}

