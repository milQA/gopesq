#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../include/pesq.h"
#include "../include/pesqpar.h"
#include "../include/dsp.h"


#define        CRITERIUM_FOR_SILENCE_OF_5_SAMPLES        500.

float Sl, Sp;



int *nr_of_hz_bands_per_bark_band;
double *centre_of_band_bark;
double *centre_of_band_hz;
double *width_of_band_bark;
double *width_of_band_hz;
double *pow_dens_correction_factor;
double *abs_thresh_power;


void make_stereo_file (char *stereo_path_name, SIGNAL_INFO *ref_info, SIGNAL_INFO *deg_info) {
    make_stereo_file2 (stereo_path_name, ref_info, deg_info-> data);
}

void make_stereo_file2 (char *stereo_path_name, SIGNAL_INFO *ref_info, float *deg) {

    long            i;
    long            h;
    short          *buffer;
    FILE           *outputFile;
    long            n;

    n = ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000) - 2 * SEARCHBUFFER * Downsample;     

    buffer = (short *) safe_malloc (2 * n * sizeof (short));
    
    if ((outputFile = fopen (stereo_path_name, "wb")) == NULL) {
        printf ("MakeStereoFile : cannot open output file %s!", stereo_path_name);
        return;
    }

    for (i = 0; i < n; i++) {
        h = (int) ref_info-> data [SEARCHBUFFER * Downsample + i] / 2;
        if (h < -32767) h = -32767;
        if (h > 32767)  h = 32767;
        h = (short) h;
        buffer [2*i] = (short) h;    
        h = (int) deg [SEARCHBUFFER * Downsample + i] / 2;
        if (h < -32767) h = -32767;
        if (h > 32767)  h = 32767;
        h = (short) h;
        buffer [2*i + 1] = (short) h;                
    }

    fwrite (buffer, sizeof (short) * 2, n, outputFile);
    
    fclose (outputFile);
    safe_free (buffer);
}

extern float InIIR_Hsos_16k [];
extern float InIIR_Hsos_8k [];
extern long InIIR_Nsos;

void select_rate( long sample_rate, long * Error_Flag, char ** Error_Type )
{
    if( Fs == sample_rate )
        return;
    if( Fs_16k == sample_rate )
    {
        Fs = Fs_16k;
        Downsample = Downsample_16k;
        InIIR_Hsos = InIIR_Hsos_16k;
        InIIR_Nsos = InIIR_Nsos_16k;
        Align_Nfft = Align_Nfft_16k;
        return;
    }
    if( Fs_8k == sample_rate )
    {
        Fs = Fs_8k;
        Downsample = Downsample_8k;
        InIIR_Hsos = InIIR_Hsos_8k;
        InIIR_Nsos = InIIR_Nsos_8k;
        Align_Nfft = Align_Nfft_8k;
        return;
    }

    (*Error_Flag) = -1;
    (*Error_Type) = "Invalid sample rate specified";    
}

int file_exist( char * fname )
{
    FILE * fp = fopen( fname, "rb" );
    if( fp == NULL )
        return 0;
    else
    {
        fclose( fp );
        return 1;
    }
}

void load_src( long * Error_Flag, char ** Error_Type,
         SIGNAL_INFO * sinfo)
{
    long name_len;
    long file_size;
    long header_size = 0;
    long Nsamples;
    long to_read;
    long read_count;
    long count;

    float *read_ptr;
    short *input_data;
    short *p_input;
    char s;
    char *p_byte;
    FILE *Src_file = fopen( sinfo-> path_name, "rb" );

    input_data = (short *) safe_malloc( 16384 * sizeof(short) );
    if( input_data == NULL )
    {
        *Error_Flag = 1;
        *Error_Type = "Could not allocate storage for file reading";
        printf ("%s!\n", *Error_Type);
        fclose( Src_file );
        return;
    }

    if( Src_file == NULL )
    {
        *Error_Flag = 1;
        *Error_Type = "Could not open source file";
        printf ("%s!\n", *Error_Type);
        safe_free( input_data );
        return;
    }

    if( fseek( Src_file, 0L, SEEK_END ) != 0 )
    {
        *Error_Flag = 1;
        *Error_Type = "Could not reach end of source file";
        safe_free( input_data );
        printf ("%s!\n", *Error_Type);
        fclose( Src_file );
        return;
    }
    file_size = ftell( Src_file );
    if( file_size < 0L )
    {
        *Error_Flag = 1;
        *Error_Type = "Could not measure length of source file";
        safe_free( input_data );
        printf ("%s!\n", *Error_Type);
        fclose( Src_file );
        return;
    }
    if( fseek( Src_file, 0L, SEEK_SET ) != 0 )
    {
        *Error_Flag = 1;
        *Error_Type = "Could not reach start of source file";
        safe_free( input_data );
        printf ("%s!\n", *Error_Type);
        fclose( Src_file );
        return;
    }
    name_len = strlen( sinfo-> path_name );
    if( name_len > 4 )
    {
        if( strcmp( sinfo-> path_name + name_len - 4, ".wav" ) == 0 )
            header_size = 22;
        if( strcmp( sinfo-> path_name + name_len - 4, ".WAV" ) == 0 )
            header_size = 22;
        if( strcmp( sinfo-> path_name + name_len - 4, ".raw" ) == 0 )
            header_size = 0;
        if( strcmp( sinfo-> path_name + name_len - 4, ".src" ) == 0 )
            header_size = 0;
    }
    if( name_len > 2 )
    {
        if( strcmp( sinfo-> path_name + name_len - 2, ".s" ) == 0 )
            header_size = 0;
    }

    /*if( header_size > 0 )
        fread( input_data, 2, header_size, Src_file );
*/
    Nsamples = (file_size / 2) - header_size;
    sinfo-> Nsamples = Nsamples + 2 * SEARCHBUFFER * Downsample;

    sinfo-> data =
        (float *) safe_malloc( (sinfo-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof(float) );
    if( sinfo-> data == NULL )
    {
        *Error_Flag = 1;
        *Error_Type = "Failed to allocate memory for source file";
        safe_free( input_data );
        printf ("%s!\n", *Error_Type);
        fclose( Src_file );
        return;
    }

    read_ptr = sinfo-> data;
    for( read_count = SEARCHBUFFER*Downsample; read_count > 0; read_count-- )
      *(read_ptr++) = 0.0f;

    to_read = Nsamples;
    while( to_read > 16384 )
    {
        read_count = fread( input_data, sizeof(short), 16384, Src_file );
        if( read_count < 16384 )
        {
            *Error_Flag = 1;
            *Error_Type = "Error reading source file.";
            printf ("%s!\n", *Error_Type);
            safe_free( input_data );
            safe_free( sinfo-> data );
            sinfo-> data = NULL;
            fclose( Src_file );
            return;
        }
        if( sinfo-> apply_swap )
        {
            p_byte = (char *)input_data;
            for( count = 0L; count < read_count; count++ )
            {
                s = p_byte[count << 1];
                p_byte[count << 1] = p_byte[(count << 1)+1];
                p_byte[(count << 1)+1] = s;
            }
        }
        to_read -= read_count;
        p_input = input_data;
        while( read_count > 0 )
        {
            read_count--;
            *(read_ptr++) = (float)(*(p_input++));
        }
    }
    read_count = fread( input_data, sizeof(short), to_read, Src_file );
    if( read_count < to_read )
    {
        *Error_Flag = 1;
        *Error_Type = "Error reading source file";
        printf ("%s!\n", *Error_Type);
        safe_free( input_data );
        safe_free( sinfo-> data );
        sinfo-> data = NULL;
        fclose( Src_file );
        return;
    }
    if( sinfo-> apply_swap )
    {
        p_byte = (char *)input_data;
        for( count = 0L; count < read_count; count++ )
        {
            s = p_byte[count << 1];
            p_byte[count << 1] = p_byte[(count << 1)+1];
            p_byte[(count << 1)+1] = s;
        }
    }
    p_input = input_data;
    while( read_count > 0 )
    {
        read_count--;
        *(read_ptr++) = (float)(*(p_input++));
    }

    for( read_count = DATAPADDING_MSECS  * (Fs / 1000) + SEARCHBUFFER * Downsample;
         read_count > 0; read_count-- )
      *(read_ptr++) = 0.0f;

    fclose( Src_file );
    safe_free( input_data );

    sinfo-> VAD = safe_malloc( sinfo-> Nsamples * sizeof(float) / Downsample );
    sinfo-> logVAD = safe_malloc( sinfo-> Nsamples * sizeof(float) / Downsample );
    if( (sinfo-> VAD == NULL) || (sinfo-> logVAD == NULL))
    {
        *Error_Flag = 1;
        *Error_Type = "Failed to allocate memory for VAD";
        printf ("%s!\n", *Error_Type);
        return;
    }
}

void alloc_other( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info, 
        long * Error_Flag, char ** Error_Type, float ** ftmp)
{
    *ftmp = (float *)safe_malloc(
       max( max(
            (*ref_info).Nsamples + DATAPADDING_MSECS  * (Fs / 1000),
            (*deg_info).Nsamples + DATAPADDING_MSECS  * (Fs / 1000) ),
           12 * Align_Nfft) * sizeof(float) );
    if( (*ftmp) == NULL )
    {
        *Error_Flag = 2;
        *Error_Type = "Failed to allocate memory for temporary storage.";
        printf ("%s!\n", *Error_Type);
        return;
    }
}

/* END OF FILE */

void input_filter(
    SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info, float * ftmp )
{
    DC_block( (*ref_info).data, (*ref_info).Nsamples );
    DC_block( (*deg_info).data, (*deg_info).Nsamples );

    apply_filters( (*ref_info).data, (*ref_info).Nsamples );
    apply_filters( (*deg_info).data, (*deg_info).Nsamples );
}

void calc_VAD( SIGNAL_INFO * sinfo )
{
    apply_VAD( sinfo, sinfo-> data, sinfo-> VAD, sinfo-> logVAD );
}

int id_searchwindows( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info )
{
    long  Utt_num = 0;
    long  count, VAD_length;
    long  this_start;
    int   speech_flag = 0;
    float VAD_value;
    long  del_deg_start;
    long  del_deg_end;

    VAD_length = ref_info-> Nsamples / Downsample;

    del_deg_start = MINUTTLENGTH - err_info-> Crude_DelayEst / Downsample;
    del_deg_end =
        ((*deg_info).Nsamples - err_info-> Crude_DelayEst) / Downsample -
        MINUTTLENGTH;

    for (count = 0; count < VAD_length; count++)
    {
        VAD_value = ref_info-> VAD [count];

        if( (VAD_value > 0.0f) && (speech_flag == 0) ) 
        {
            speech_flag = 1;
            this_start = count;
            err_info-> UttSearch_Start [Utt_num] = count - SEARCHBUFFER;
            if( err_info-> UttSearch_Start [Utt_num] < 0 )
                err_info-> UttSearch_Start [Utt_num] = 0;
        }

        if( ((VAD_value == 0.0f) || (count == (VAD_length-1))) &&
            (speech_flag == 1) ) 
        {
            speech_flag = 0;
            err_info-> UttSearch_End [Utt_num] = count + SEARCHBUFFER;
            if( err_info-> UttSearch_End [Utt_num] > VAD_length - 1 )
                err_info-> UttSearch_End [Utt_num] = VAD_length -1;

            if( ((count - this_start) >= MINUTTLENGTH) &&
                (this_start < del_deg_end) &&
                (count > del_deg_start) )
                Utt_num++;            
        }
    }

    err_info-> Nutterances = Utt_num;
    return Utt_num;
} 

void id_utterances( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info )
{
    long  Utt_num = 0;
    long  Largest_uttsize = 0;
    long  count, VAD_length;
    int   speech_flag = 0;
    float VAD_value;
    long  this_start;
    long  last_end;
    long  del_deg_start;
    long  del_deg_end;

    VAD_length = ref_info-> Nsamples / Downsample;

    del_deg_start = MINUTTLENGTH - err_info-> Crude_DelayEst / Downsample;
    del_deg_end =
        ((*deg_info).Nsamples - err_info-> Crude_DelayEst) / Downsample -
        MINUTTLENGTH;

    for (count = 0; count < VAD_length ; count++)
    {
        VAD_value = ref_info-> VAD [count];
        if( (VAD_value > 0.0f) && (speech_flag == 0) ) 
        {
            speech_flag = 1;
            this_start = count;
            err_info-> Utt_Start [Utt_num] = count;
        }

        if( ((VAD_value == 0.0f) || (count == (VAD_length-1))) &&
            (speech_flag == 1) ) 
        {
            speech_flag = 0;
            err_info-> Utt_End [Utt_num] = count;

            if( ((count - this_start) >= MINUTTLENGTH) &&
                (this_start < del_deg_end) &&
                (count > del_deg_start) )
                Utt_num++;            
        }
    }

    err_info-> Utt_Start [0] = SEARCHBUFFER;
    err_info-> Utt_End [err_info-> Nutterances-1] = (VAD_length - SEARCHBUFFER);
    
    for (Utt_num = 1; Utt_num < err_info-> Nutterances; Utt_num++ )
    {
        this_start = err_info-> Utt_Start [Utt_num];
        last_end = err_info-> Utt_End [Utt_num - 1];
        count = (this_start + last_end) / 2;
        err_info-> Utt_Start [Utt_num] = count;
        err_info-> Utt_End [Utt_num - 1] = count;
    }

    this_start = (err_info-> Utt_Start [0] * Downsample) + err_info-> Utt_Delay [0];
    if( this_start < (SEARCHBUFFER * Downsample) )
    {
        count = SEARCHBUFFER +
                (Downsample - 1 - err_info-> Utt_Delay [0]) / Downsample;
        err_info-> Utt_Start [0] = count;
    }
    last_end = (err_info-> Utt_End [err_info-> Nutterances-1] * Downsample) +
               err_info-> Utt_Delay [err_info-> Nutterances-1];
    if( last_end > ((*deg_info).Nsamples - SEARCHBUFFER * Downsample) )
    {
        count = ( (*deg_info).Nsamples -
                  err_info-> Utt_Delay [err_info-> Nutterances-1] ) / Downsample -
                SEARCHBUFFER;
        err_info-> Utt_End [err_info-> Nutterances-1] = count;
    }

    for (Utt_num = 1; Utt_num < err_info-> Nutterances; Utt_num++ )
    {
        this_start =
            (err_info-> Utt_Start [Utt_num] * Downsample) +
            err_info-> Utt_Delay [Utt_num];
        last_end =
            (err_info-> Utt_End [Utt_num - 1] * Downsample) +
            err_info-> Utt_Delay [Utt_num - 1];
        if( this_start < last_end )
        {
            count = (this_start + last_end) / 2;
            this_start =
                (Downsample - 1 + count - err_info-> Utt_Delay [Utt_num]) / Downsample;
            last_end =
               (count - err_info-> Utt_Delay [Utt_num - 1]) / Downsample;
            err_info-> Utt_Start [Utt_num] = this_start;
            err_info-> Utt_End [Utt_num - 1] = last_end;
        }
    }

    for (Utt_num = 0; Utt_num < err_info-> Nutterances; Utt_num++ )
        if( (err_info-> Utt_End [Utt_num] - err_info-> Utt_Start [Utt_num])
             > Largest_uttsize )
            Largest_uttsize = 
                err_info-> Utt_End [Utt_num] - err_info-> Utt_Start [Utt_num];

    err_info-> Largest_uttsize = Largest_uttsize;
}

void utterance_split( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info, float * ftmp )
{
    long Utt_id;
    long Utt_DelayEst;
    long Utt_Delay;
    float Utt_DelayConf;
    long Utt_Start;
    long Utt_End;
    long Utt_SpeechStart;
    long Utt_SpeechEnd;
    long Utt_Len;
    long step;
    long Best_ED1, Best_ED2;
    long Best_D1, Best_D2;
    float Best_DC1, Best_DC2;
    long Best_BP;
    long Largest_uttsize = 0;

    Utt_id = 0;
    while( (Utt_id < err_info-> Nutterances) &&
           (err_info-> Nutterances < MAXNUTTERANCES) )
    {
        Utt_DelayEst = err_info-> Utt_DelayEst [Utt_id];
        Utt_Delay = err_info-> Utt_Delay [Utt_id];
        Utt_DelayConf = err_info-> Utt_DelayConf [Utt_id];
        Utt_Start = err_info-> Utt_Start [Utt_id];
        Utt_End = err_info-> Utt_End [Utt_id];

        Utt_SpeechStart = Utt_Start;
        while( (Utt_SpeechStart < Utt_End) && (ref_info-> VAD [Utt_SpeechStart] <= 0.0f) )
            Utt_SpeechStart++;
        Utt_SpeechEnd = Utt_End;
        while( (Utt_SpeechEnd > Utt_Start) && (ref_info-> VAD [Utt_SpeechEnd] <= 0.0f) )
            Utt_SpeechEnd--;
        Utt_SpeechEnd++;
        Utt_Len = Utt_SpeechEnd - Utt_SpeechStart;

        if( Utt_Len >= 200 )
        {
            split_align( ref_info, deg_info, err_info, ftmp,
                Utt_Start, Utt_SpeechStart, Utt_SpeechEnd, Utt_End,
                Utt_DelayEst, Utt_DelayConf,
                &Best_ED1, &Best_D1, &Best_DC1,
                &Best_ED2, &Best_D2, &Best_DC2,
                &Best_BP );

            if( (Best_DC1 > Utt_DelayConf) && (Best_DC2 > Utt_DelayConf) )
            {
                for (step = err_info-> Nutterances-1; step > Utt_id; step-- )
                {
                    err_info-> Utt_DelayEst [step +1] = err_info-> Utt_DelayEst [step];
                    err_info-> Utt_Delay [step +1] = err_info-> Utt_Delay [step];
                    err_info-> Utt_DelayConf [step +1] = err_info-> Utt_DelayConf [step];
                    err_info-> Utt_Start [step +1] = err_info-> Utt_Start [step];
                    err_info-> Utt_End [step +1] = err_info-> Utt_End [step];
                    err_info-> UttSearch_Start [step +1] = err_info-> Utt_Start [step];
                    err_info-> UttSearch_End [step +1] = err_info-> Utt_End [step];
                }
                err_info-> Nutterances++;

                err_info-> Utt_DelayEst [Utt_id] = Best_ED1;
                err_info-> Utt_Delay [Utt_id] = Best_D1;
                err_info-> Utt_DelayConf [Utt_id] = Best_DC1;

                err_info-> Utt_DelayEst [Utt_id +1] = Best_ED2;
                err_info-> Utt_Delay [Utt_id +1] = Best_D2;
                err_info-> Utt_DelayConf [Utt_id +1] = Best_DC2;

                err_info-> UttSearch_Start [Utt_id +1] = err_info-> UttSearch_Start [Utt_id];
                err_info-> UttSearch_End [Utt_id +1] = err_info-> UttSearch_End [Utt_id];

                if( Best_D2 < Best_D1 )
                {
                    err_info-> Utt_Start [Utt_id] = Utt_Start;
                    err_info-> Utt_End [Utt_id] = Best_BP;
                    err_info-> Utt_Start [Utt_id +1] = Best_BP;
                    err_info-> Utt_End [Utt_id +1] = Utt_End;
                }
                else
                {
                    err_info-> Utt_Start [Utt_id] = Utt_Start;
                    err_info-> Utt_End [Utt_id] = Best_BP + (Best_D2 - Best_D1) / (2 * Downsample);
                    err_info-> Utt_Start [Utt_id +1] = Best_BP - (Best_D2 - Best_D1) / (2 * Downsample);
                    err_info-> Utt_End [Utt_id +1] = Utt_End;
                }

                if( (err_info-> Utt_Start [Utt_id] - SEARCHBUFFER) * Downsample + Best_D1 < 0 )
                    err_info-> Utt_Start [Utt_id] =
                        SEARCHBUFFER + (Downsample - 1 - Best_D1) / Downsample;

                if( (err_info-> Utt_End [Utt_id +1] * Downsample + Best_D2) >
                    ((*deg_info).Nsamples - SEARCHBUFFER * Downsample) )
                    err_info-> Utt_End [Utt_id +1] =
                        ((*deg_info).Nsamples - Best_D2) / Downsample - SEARCHBUFFER;

            }
            else Utt_id++;
        }
        else Utt_id++;
    }

    for (Utt_id = 0; Utt_id < err_info-> Nutterances; Utt_id++ )
        if( (err_info-> Utt_End [Utt_id] - err_info-> Utt_Start [Utt_id])
             > Largest_uttsize )
            Largest_uttsize = 
                err_info-> Utt_End [Utt_id] - err_info-> Utt_Start [Utt_id];

    err_info-> Largest_uttsize = Largest_uttsize;
}

void utterance_locate( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info, float * ftmp )
{    
    long Utt_id;
    
    id_searchwindows( ref_info, deg_info, err_info );

    for (Utt_id = 0; Utt_id < err_info-> Nutterances; Utt_id++)
    {
        crude_align( ref_info, deg_info, err_info, Utt_id, ftmp);
        time_align(ref_info, deg_info, err_info, Utt_id, ftmp );
    }

    id_utterances( ref_info, deg_info, err_info );

    utterance_split( ref_info, deg_info, err_info, ftmp );   
}


void short_term_fft (int Nf, SIGNAL_INFO *info, float *window, long start_sample, float *hz_spectrum, float *fft_tmp) {
    int n, k;        

    for (n = 0; n < Nf; n++ )
    {
        fft_tmp [n] = info-> data [start_sample + n] * window [n];
    }
    RealFFT(fft_tmp, Nf);

    for (k = 0; k < Nf / 2; k++ ) 
    {
        hz_spectrum [k] = fft_tmp [k << 1] * fft_tmp [k << 1] + fft_tmp [1 + (k << 1)] * fft_tmp [1 + (k << 1)];
    }    

    hz_spectrum [0] = 0;
}

void freq_warping (int number_of_hz_bands, float *hz_spectrum, int Nb, float *pitch_pow_dens, long frame) {

    int        hz_band = 0;
    int        bark_band;
    double    sum;

    for (bark_band = 0; bark_band < Nb; bark_band++) {
        int n = nr_of_hz_bands_per_bark_band [bark_band];
        int i;

        sum = 0;
        for (i = 0; i < n; i++) {
            sum += hz_spectrum [hz_band++];
        }
        
        sum *= pow_dens_correction_factor [bark_band];
        sum *= Sp;
        pitch_pow_dens [frame * Nb + bark_band] = (float) sum;
    }
}

float total_audible (int frame, float *pitch_pow_dens, float factor) {
    int        band;
    float     h, threshold;
    double  result;
    
    result = 0.;
    for (band= 1; band< Nb; band++) {
        h = pitch_pow_dens [frame * Nb + band];
        threshold = (float) (factor * abs_thresh_power [band]);
        if (h > threshold) {
            result += h;
        }
    }
    return (float) result;
}

void time_avg_audible_of (int number_of_frames, int *silent, float *pitch_pow_dens, float *avg_pitch_pow_dens, int total_number_of_frames) 
{
    int    frame;
    int    band;

    for (band = 0; band < Nb; band++) {
        double result = 0;
        for (frame = 0; frame < number_of_frames; frame++) {
            if (!silent [frame]) {
                float h = pitch_pow_dens [frame * Nb + band];
                if (h > 100 * abs_thresh_power [band]) {
                    result += h;
                }
            }
        }

        avg_pitch_pow_dens [band] = (float) (result / total_number_of_frames);
    }
}            

void freq_resp_compensation (int number_of_frames, float *pitch_pow_dens_ref, float *avg_pitch_pow_dens_ref, float *avg_pitch_pow_dens_deg, float constant)
{
    int band;

    for (band = 0; band < Nb; band++) {
        float    x = (avg_pitch_pow_dens_deg [band] + constant) / (avg_pitch_pow_dens_ref [band] + constant);
        int        frame;

        if (x > (float) 100.0) {x = (float) 100.0;} 
        if (x < (float) 0.01) {x = (float) 0.01;}   

        for (frame = 0; frame < number_of_frames; frame++) {        
            pitch_pow_dens_ref [frame * Nb + band] *= x;
        }        
    }
}

#define ZWICKER_POWER       0.23 

void intensity_warping_of (float *loudness_dens, int frame, float *pitch_pow_dens)
{
    int        band;
    float    h;
    double    modified_zwicker_power;

    for (band = 0; band < Nb; band++) {
        float threshold = (float) abs_thresh_power [band];
        float input = pitch_pow_dens [frame * Nb + band];

        if (centre_of_band_bark [band] < (float) 4) {
            h =  (float) 6 / ((float) centre_of_band_bark [band] + (float) 2);
        } else {
            h = (float) 1;
        }
        if (h > (float) 2) {h = (float) 2;}
        h = (float) pow (h, (float) 0.15); 
        modified_zwicker_power = ZWICKER_POWER * h;

        if (input > threshold) {
            loudness_dens [band] = (float) (pow (threshold / 0.5, modified_zwicker_power)
                                                    * (pow (0.5 + 0.5 * input / threshold, modified_zwicker_power) - 1));
        } else {
            loudness_dens [band] = 0;
        }

        loudness_dens [band] *= (float) Sl;
    }    
}

float pseudo_Lp (int n, float *x, float p) {   
    double totalWeight = 0;
    double result = 0;
    int    band;

    for (band = 1; band < Nb; band++) {
        float h = (float) fabs (x [band]);        
        float w = (float) width_of_band_bark [band];
        float prod = h * w;

        result += pow (prod, p);
        totalWeight += w;
    }

    result /= totalWeight;
    result = pow (result, 1/p);
    result *= totalWeight;
    
    return (float) result;
}  
void multiply_with_asymmetry_factor (float      *disturbance_dens, 
                                     int         frame, 
                                     const float   * const pitch_pow_dens_ref, 
                                     const float   * const pitch_pow_dens_deg) 
{
    int   i;
    float ratio, h;

    for (i = 0; i < Nb; i++) {
        ratio = (pitch_pow_dens_deg [frame * Nb + i] + (float) 50)
                  / (pitch_pow_dens_ref [frame * Nb + i] + (float) 50);

        h = (float) pow (ratio, (float) 1.2);    
        if (h > (float) 12) {h = (float) 12;}
        if (h < (float) 3) {h = (float) 0.0;}

        disturbance_dens [i] *= h;
    }
}

double pow_of (const float * const x, long start_sample, long stop_sample, long divisor) {
    long    i;
    double  power = 0;

    if (start_sample < 0) {
        exit (1);
    }

    if (start_sample > stop_sample) {
        exit (1);
    }

    for (i = start_sample; i < stop_sample; i++) {
        float h = x [i];
        power += h * h;        
    }
    
    power /= divisor;
    return power;
}


int compute_delay (long              start_sample, 
                   long                 stop_sample, 
                   long                 search_range, 
                   float            *time_series1, 
                   float            *time_series2,
                   float            *max_correlation) {

    double            power1, power2, normalization;
    long            i;
    float           *x1, *x2, *y;
    double            h;
    long            n = stop_sample - start_sample;   
    long            power_of_2 = nextpow2 (2 * n);
    long            best_delay;

    power1 = pow_of (time_series1, start_sample, stop_sample, stop_sample - start_sample) * (double) n/(double) power_of_2;
    power2 = pow_of (time_series2, start_sample, stop_sample, stop_sample - start_sample) * (double) n/(double) power_of_2;
    normalization = sqrt (power1 * power2);

    if ((power1 <= 1E-6) || (power2 <= 1E-6)) {
        *max_correlation = 0;
        return 0;
    }

    x1 = (float *) safe_malloc ((power_of_2 + 2) * sizeof (float));;
    x2 = (float *) safe_malloc ((power_of_2 + 2) * sizeof (float));;
    y = (float *) safe_malloc ((power_of_2 + 2) * sizeof (float));;
    
    for (i = 0; i < power_of_2 + 2; i++) {
        x1 [i] = 0.;
        x2 [i] = 0.;
        y [i] = 0.;
    }

    for (i = 0; i < n; i++) {
        x1 [i] = (float) fabs (time_series1 [i + start_sample]);
        x2 [i] = (float) fabs (time_series2 [i + start_sample]);
    }

    RealFFT (x1, power_of_2);
    RealFFT (x2, power_of_2);

    for (i = 0; i <= power_of_2 / 2; i++) { 
        x1 [2 * i] /= power_of_2;
        x1 [2 * i + 1] /= power_of_2;                
    }

    for (i = 0; i <= power_of_2 / 2; i++) { 
        y [2*i] = x1 [2*i] * x2 [2*i] + x1 [2*i + 1] * x2 [2*i + 1];
        y [2*i + 1] = -x1 [2*i + 1] * x2 [2*i] + x1 [2*i] * x2 [2*i + 1];
    }    
  
    RealIFFT (y, power_of_2);

    best_delay = 0;
    *max_correlation = 0;

    for (i = -search_range; i <= -1; i++) {
        h = (float) fabs (y [(i + power_of_2)]) / normalization;
        if (fabs (h) > (double) *max_correlation) {
            *max_correlation = (float) fabs (h);
            best_delay= i;
        }
    }

    for (i = 0; i < search_range; i++) {
        h = (float) fabs (y [i]) / normalization;
        if (fabs (h) > (double) *max_correlation) {
            *max_correlation = (float) fabs (h);
            best_delay= i;
        }
    }

    safe_free (x1);
    safe_free (x2);
    safe_free (y);
    
    return best_delay;
}

    
#define NUMBER_OF_PSQM_FRAMES_PER_SYLLABE       20
 

float Lpq_weight (int         start_frame,
                  int         stop_frame,
                  float         power_syllable,
                  float      power_time,
                  float        *frame_disturbance,
                  float        *time_weight) {

    double    result_time= 0;
    double  total_time_weight_time = 0;
    int        start_frame_of_syllable;
    
    for (start_frame_of_syllable = start_frame; 
         start_frame_of_syllable <= stop_frame; 
         start_frame_of_syllable += NUMBER_OF_PSQM_FRAMES_PER_SYLLABE/2) {

        double  result_syllable = 0;
        int     count_syllable = 0;
        int     frame;

        for (frame = start_frame_of_syllable;
             frame < start_frame_of_syllable + NUMBER_OF_PSQM_FRAMES_PER_SYLLABE;
             frame++) {
            if (frame <= stop_frame) {
                float h = frame_disturbance [frame];
                result_syllable +=  pow (h, power_syllable); 
            }
            count_syllable++;                
        }

        result_syllable /= count_syllable;
        result_syllable = pow (result_syllable, (double) 1/power_syllable);        
    
        result_time+=  pow (time_weight [start_frame_of_syllable - start_frame] * result_syllable, power_time); 
        total_time_weight_time += pow (time_weight [start_frame_of_syllable - start_frame], power_time);
    }

    result_time /= total_time_weight_time;
    result_time= pow (result_time, (float) 1 / power_time);

    return (float) result_time;
}

void set_to_sine (SIGNAL_INFO *info, float amplitude, float omega) {
    long i;

    for (i = 0; i < info-> Nsamples; i++) {
        info-> data [i] = amplitude * (float) sin (omega * i);
    }
}

float maximum_of (float *x, long start, long stop) {
    long i;
    float result = -1E20f;

    for (i = start; i < stop; i++) {
        if (result < x [i]) {
            result = x [i];
        }
    }

    return result;
}

float integral_of (float *x, long frames_after_start) {
    double result = 0;
    int    band;

    for (band = 1; band < Nb; band++) {
        result += x [frames_after_start * Nb + band] * width_of_band_bark [band];        
    }
    return (float) result;    


    return (float) result;
}

#define DEBUG_FR    0

void pesq_psychoacoustic_model(SIGNAL_INFO    * ref_info, 
                                 SIGNAL_INFO    * deg_info,
                               ERROR_INFO    * err_info, 
                               float        * ftmp)
{

    long    maxNsamples = max (ref_info-> Nsamples, deg_info-> Nsamples);
    long    Nf = Downsample * 8L;
    long    start_frame, stop_frame;
    long    samples_to_skip_at_start, samples_to_skip_at_end;
    float   sum_of_5_samples;
    long    n, i;
    float   power_ref, power_deg;
    long    frame;
    float   *fft_tmp;
    float    *hz_spectrum_ref, *hz_spectrum_deg;
    float   *pitch_pow_dens_ref, *pitch_pow_dens_deg;
    float    *loudness_dens_ref, *loudness_dens_deg;
    float   *avg_pitch_pow_dens_ref, *avg_pitch_pow_dens_deg;
    float    *deadzone;
    float   *disturbance_dens, *disturbance_dens_asym_add;
    float     total_audible_pow_ref, total_audible_pow_deg;
    int        *silent;
    float    oldScale, scale;
    int     *frame_was_skipped;
    float   *frame_disturbance;
    float   *frame_disturbance_asym_add;
    float   *total_power_ref;
    int         utt;
    
#ifdef CALIBRATE
    int     periodInSamples;
    int     numberOfPeriodsPerFrame;
    float   omega; 
#endif

    float   peak;

#define    MAX_NUMBER_OF_BAD_INTERVALS        1000

    int        *frame_is_bad; 
    int        *smeared_frame_is_bad; 
    int         start_frame_of_bad_interval [MAX_NUMBER_OF_BAD_INTERVALS];    
    int         stop_frame_of_bad_interval [MAX_NUMBER_OF_BAD_INTERVALS];    
    int         start_sample_of_bad_interval [MAX_NUMBER_OF_BAD_INTERVALS];    
    int         stop_sample_of_bad_interval [MAX_NUMBER_OF_BAD_INTERVALS];   
    int         number_of_samples_in_bad_interval [MAX_NUMBER_OF_BAD_INTERVALS];    
    int         delay_in_samples_in_bad_interval  [MAX_NUMBER_OF_BAD_INTERVALS];    
    int         number_of_bad_intervals= 0;
    int         search_range_in_samples;
    int         bad_interval;
    float *untweaked_deg = NULL;
    float *tweaked_deg = NULL;
    float *doubly_tweaked_deg = NULL;
    int         there_is_a_bad_frame = FALSE;
    float    *time_weight;
    float    d_indicator, a_indicator;
    int      nn;

    float Whanning [Nfmax];

    for (n = 0L; n < Nf; n++ ) {
        Whanning [n] = (float)(0.5 * (1.0 - cos((TWOPI * n) / Nf)));
    }

    switch (Fs) {
    case 8000:
        Nb = 42;
        Sl = (float) Sl_8k;
        Sp = (float) Sp_8k;
        nr_of_hz_bands_per_bark_band = nr_of_hz_bands_per_bark_band_8k;
        centre_of_band_bark = centre_of_band_bark_8k;
        centre_of_band_hz = centre_of_band_hz_8k;
        width_of_band_bark = width_of_band_bark_8k;
        width_of_band_hz = width_of_band_hz_8k;
        pow_dens_correction_factor = pow_dens_correction_factor_8k;
        abs_thresh_power = abs_thresh_power_8k;
        break;
    case 16000:
        Nb = 49;
        Sl = (float) Sl_16k;
        Sp = (float) Sp_16k;
        nr_of_hz_bands_per_bark_band = nr_of_hz_bands_per_bark_band_16k;
        centre_of_band_bark = centre_of_band_bark_16k;
        centre_of_band_hz = centre_of_band_hz_16k;
        width_of_band_bark = width_of_band_bark_16k;
        width_of_band_hz = width_of_band_hz_16k;
        pow_dens_correction_factor = pow_dens_correction_factor_16k;
        abs_thresh_power = abs_thresh_power_16k;
        break;
    default:
        printf ("Invalid sample frequency!\n");
        exit (1);
    }

    samples_to_skip_at_start = 0;
    do {
        sum_of_5_samples= (float) 0;
        for (i = 0; i < 5; i++) {
            sum_of_5_samples += (float) fabs (ref_info-> data [SEARCHBUFFER * Downsample + samples_to_skip_at_start + i]);
        }
        if (sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) {
            samples_to_skip_at_start++;         
        }        
    } while ((sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) 
            && (samples_to_skip_at_start < maxNsamples / 2));
    
    samples_to_skip_at_end = 0;
    do {
        sum_of_5_samples= (float) 0;
        for (i = 0; i < 5; i++) {
            sum_of_5_samples += (float) fabs (ref_info-> data [maxNsamples - SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000) - 1 - samples_to_skip_at_end - i]);
        }
        if (sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) {
            samples_to_skip_at_end++;         
        }        
    } while ((sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) 
        && (samples_to_skip_at_end < maxNsamples / 2));
       
    start_frame = samples_to_skip_at_start / (Nf /2);
    stop_frame = (maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000) - samples_to_skip_at_end) / (Nf /2) - 1; 

    power_ref = (float) pow_of (ref_info-> data, 
                                SEARCHBUFFER * Downsample, 
                                maxNsamples - SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000),
                                maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000)); 
    power_deg = (float) pow_of (deg_info-> data, 
                                SEARCHBUFFER * Downsample, 
                                maxNsamples - SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000),
                                maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000));

    fft_tmp                = (float *) safe_malloc ((Nf + 2) * sizeof (float));
    hz_spectrum_ref        = (float *) safe_malloc ((Nf / 2) * sizeof (float));
    hz_spectrum_deg        = (float *) safe_malloc ((Nf / 2) * sizeof (float));
    
    frame_is_bad        = (int *) safe_malloc ((stop_frame + 1) * sizeof (int)); 
    smeared_frame_is_bad=(int *) safe_malloc ((stop_frame + 1) * sizeof (int)); 
    
    silent                = (int *) safe_malloc ((stop_frame + 1) * sizeof (int));
    
    pitch_pow_dens_ref    = (float *) safe_malloc ((stop_frame + 1) * Nb * sizeof (float));
    pitch_pow_dens_deg    = (float *) safe_malloc ((stop_frame + 1) * Nb * sizeof (float));
    
    frame_was_skipped    = (int *) safe_malloc ((stop_frame + 1) * sizeof (int));

    frame_disturbance    = (float *) safe_malloc ((stop_frame + 1) * sizeof (float));
    frame_disturbance_asym_add    = (float *) safe_malloc ((stop_frame + 1) * sizeof (float));

    avg_pitch_pow_dens_ref = (float *) safe_malloc (Nb * sizeof (float));
    avg_pitch_pow_dens_deg = (float *) safe_malloc (Nb * sizeof (float));
    loudness_dens_ref    = (float *) safe_malloc (Nb * sizeof (float));
    loudness_dens_deg    = (float *) safe_malloc (Nb * sizeof (float));;
    deadzone                = (float *) safe_malloc (Nb * sizeof (float));;
    disturbance_dens    = (float *) safe_malloc (Nb * sizeof (float));
    disturbance_dens_asym_add = (float *) safe_malloc (Nb * sizeof (float));    

    time_weight            = (float *) safe_malloc ((stop_frame + 1) * sizeof (float));
    total_power_ref     = (float *) safe_malloc ((stop_frame + 1) * sizeof (float));
        
#ifdef CALIBRATE
    periodInSamples = Fs / 1000;
    numberOfPeriodsPerFrame = Nf / periodInSamples;
    omega = (float) (TWOPI / periodInSamples);    
    peak;

    set_to_sine (ref_info, (float) 29.54, (float) omega);    
#endif

    for (frame = 0; frame <= stop_frame; frame++) {
        int start_sample_ref = SEARCHBUFFER * Downsample + frame * Nf / 2;
        int start_sample_deg;
        int delay;    

        short_term_fft (Nf, ref_info, Whanning, start_sample_ref, hz_spectrum_ref, fft_tmp);
        
        if (err_info-> Nutterances < 1) {
            printf ("Processing error!\n");
            exit (1);
        }

        utt = err_info-> Nutterances - 1;
        while ((utt >= 0) && (err_info-> Utt_Start [utt] * Downsample > start_sample_ref)) {
            utt--;
        }
        if (utt >= 0) {
            delay = err_info-> Utt_Delay [utt];
        } else {
            delay = err_info-> Utt_Delay [0];        
        }
        start_sample_deg = start_sample_ref + delay;         

        if ((start_sample_deg > 0) && (start_sample_deg + Nf < maxNsamples + DATAPADDING_MSECS  * (Fs / 1000))) {
            short_term_fft (Nf, deg_info, Whanning, start_sample_deg, hz_spectrum_deg, fft_tmp);            
        } else {
            for (i = 0; i < Nf / 2; i++) {
                hz_spectrum_deg [i] = 0;
            }
        }

        freq_warping (Nf / 2, hz_spectrum_ref, Nb, pitch_pow_dens_ref, frame);

        peak = maximum_of (pitch_pow_dens_ref, 0, Nb);    

        freq_warping (Nf / 2, hz_spectrum_deg, Nb, pitch_pow_dens_deg, frame);
    
        total_audible_pow_ref = total_audible (frame, pitch_pow_dens_ref, 1E2);
        total_audible_pow_deg = total_audible (frame, pitch_pow_dens_deg, 1E2);        

        silent [frame] = (total_audible_pow_ref < 1E7);     
    }    

    time_avg_audible_of (stop_frame + 1, silent, pitch_pow_dens_ref, avg_pitch_pow_dens_ref, (maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000)) / (Nf / 2) - 1);
    time_avg_audible_of (stop_frame + 1, silent, pitch_pow_dens_deg, avg_pitch_pow_dens_deg, (maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000)) / (Nf / 2) - 1);
    
#ifndef CALIBRATE
    freq_resp_compensation (stop_frame + 1, pitch_pow_dens_ref, avg_pitch_pow_dens_ref, avg_pitch_pow_dens_deg, 1000);
#endif
    
    oldScale = 1;
    for (frame = 0; frame <= stop_frame; frame++) {
        int band;

        total_audible_pow_ref = total_audible (frame, pitch_pow_dens_ref, 1);
        total_audible_pow_deg = total_audible (frame, pitch_pow_dens_deg, 1);        
        total_power_ref [frame] = total_audible_pow_ref;

        scale = (total_audible_pow_ref + (float) 5E3) / (total_audible_pow_deg + (float) 5E3);
                
        if (frame > 0) {
            scale = (float) 0.2 * oldScale + (float) 0.8*scale;
        }
        oldScale = scale;

#define MAX_SCALE   5.0
    
        if (scale > (float) MAX_SCALE) scale = (float) MAX_SCALE;

#define MIN_SCALE   3E-4
    
        if (scale < (float) MIN_SCALE) {
            scale = (float) MIN_SCALE;            
        }

        for (band = 0; band < Nb; band++) {
            pitch_pow_dens_deg [frame * Nb + band] *= scale;
        }

        intensity_warping_of (loudness_dens_ref, frame, pitch_pow_dens_ref); 
        intensity_warping_of (loudness_dens_deg, frame, pitch_pow_dens_deg); 
    
        for (band = 0; band < Nb; band++) {
            disturbance_dens [band] = loudness_dens_deg [band] - loudness_dens_ref [band];
        }
        
        for (band = 0; band < Nb; band++) {
            deadzone [band] = min (loudness_dens_deg [band], loudness_dens_ref [band]);    
            deadzone [band] *= 0.25;
        }
        
        for (band = 0; band < Nb; band++) {
            float d = disturbance_dens [band];
            float m = deadzone [band];
                
            if (d > m) {
                disturbance_dens [band] -= m;
            } else {
                if (d < -m) {
                    disturbance_dens [band] += m;
                } else {
                    disturbance_dens [band] = 0;
                }
            }
        }    

        frame_disturbance [frame] = pseudo_Lp (Nb, disturbance_dens, D_POW_F);    

#define THRESHOLD_BAD_FRAMES   30

        if (frame_disturbance [frame] > THRESHOLD_BAD_FRAMES) 
        {
            there_is_a_bad_frame = TRUE;
        }

        multiply_with_asymmetry_factor (disturbance_dens, frame, pitch_pow_dens_ref, pitch_pow_dens_deg);
    
        frame_disturbance_asym_add [frame] = pseudo_Lp (Nb, disturbance_dens, A_POW_F);    
    }

    for (frame = 0; frame <= stop_frame; frame++) {
        frame_was_skipped [frame] = FALSE;
    }

    for (utt = 1; utt < err_info-> Nutterances; utt++) {
        int frame1 = (int) floor (((err_info-> Utt_Start [utt] - SEARCHBUFFER ) * Downsample + err_info-> Utt_Delay [utt]) / (Nf / 2));
        int j = (int) floor ((err_info-> Utt_End [utt-1] - SEARCHBUFFER) * Downsample + err_info-> Utt_Delay [utt-1]) / (Nf / 2);
        int delay_jump = err_info-> Utt_Delay [utt] - err_info-> Utt_Delay [utt-1];

        if (frame1 > j) {
            frame1 = j;
        }
        
        if (frame1 < 0) {
            frame1 = 0;
        }
            
        if (delay_jump < -(int) (Nf / 2)) {

            int frame2 = (int) ((err_info-> Utt_Start [utt] - SEARCHBUFFER) * Downsample + max (0, fabs (delay_jump))) / (Nf / 2) + 1; 
            
            for (frame = frame1; frame <= frame2; frame++)  {
                if (frame < stop_frame) {
                    frame_was_skipped [frame] = TRUE;

                    frame_disturbance [frame] = 0;
                    frame_disturbance_asym_add [frame] = 0;
                }
            } 
        }    
    }
    
    nn = DATAPADDING_MSECS  * (Fs / 1000) + maxNsamples;

    tweaked_deg = (float *) safe_malloc (nn * sizeof (float));

    for (i = 0; i < nn; i++) {
        tweaked_deg [i] = 0;
    } 

    for (i = SEARCHBUFFER * Downsample; i < nn - SEARCHBUFFER * Downsample; i++) {
        int  utt = err_info-> Nutterances - 1;
        long delay, j;

        while ((utt >= 0) && (err_info-> Utt_Start [utt] * Downsample > i)) {
            utt--;
        }
        if (utt >= 0) {
            delay = err_info-> Utt_Delay [utt];
        } else {
            delay = err_info-> Utt_Delay [0];        
        }
            
        j = i + delay;
        if (j < SEARCHBUFFER * Downsample) {
            j = SEARCHBUFFER * Downsample;
        }
        if (j >= nn - SEARCHBUFFER * Downsample) {
            j = nn - SEARCHBUFFER * Downsample - 1;
        }
        tweaked_deg [i] = deg_info-> data [j];
    }

    if (there_is_a_bad_frame) {        
        
        for (frame = 0; frame <= stop_frame; frame++) 
        {  
            frame_is_bad [frame] = (frame_disturbance [frame] > THRESHOLD_BAD_FRAMES);       

            smeared_frame_is_bad [frame] = FALSE;
        }
        frame_is_bad [0] = FALSE;

    #define SMEAR_RANGE 2
        
        for (frame = SMEAR_RANGE; frame < stop_frame - SMEAR_RANGE; frame++) {    
            long max_itself_and_left = frame_is_bad [frame];
            long max_itself_and_right = frame_is_bad [frame];
            long mini, i;

            for (i = -SMEAR_RANGE; i <= 0; i++) {
                if (max_itself_and_left < frame_is_bad [frame  + i]) {
                    max_itself_and_left = frame_is_bad [frame  + i];
                }
            }
        
            for (i = 0; i <= SMEAR_RANGE; i++) {
                if (max_itself_and_right < frame_is_bad [frame + i]) {
                    max_itself_and_right = frame_is_bad [frame + i];
                }
            }

            mini = max_itself_and_left;
            if (mini > max_itself_and_right) {
                mini = max_itself_and_right;
            }

            smeared_frame_is_bad [frame] = mini;
        }
   
#define MINIMUM_NUMBER_OF_BAD_FRAMES_IN_BAD_INTERVAL    5

        number_of_bad_intervals = 0;    
        frame = 0; 
        while (frame <= stop_frame) {

            while ((frame <= stop_frame) && (!smeared_frame_is_bad [frame])) {
                frame++; 
            }

            if (frame <= stop_frame) { 
                start_frame_of_bad_interval [number_of_bad_intervals] = frame;

                while ((frame <= stop_frame) && (smeared_frame_is_bad [frame])) {
                    frame++; 
                }
            
                if (frame <= stop_frame) {
                    stop_frame_of_bad_interval [number_of_bad_intervals] = frame; 

                    if (stop_frame_of_bad_interval [number_of_bad_intervals] - start_frame_of_bad_interval [number_of_bad_intervals] >= MINIMUM_NUMBER_OF_BAD_FRAMES_IN_BAD_INTERVAL) {
                        number_of_bad_intervals++; 
                    }
                }
            }
        }

        for (bad_interval = 0; bad_interval < number_of_bad_intervals; bad_interval++) {
            start_sample_of_bad_interval [bad_interval] =  start_frame_of_bad_interval [bad_interval] * (Nf / 2) + SEARCHBUFFER * Downsample;
            stop_sample_of_bad_interval [bad_interval] =  stop_frame_of_bad_interval [bad_interval] * (Nf / 2) + Nf + SEARCHBUFFER* Downsample;
            if (stop_frame_of_bad_interval [bad_interval] > stop_frame) {
                stop_frame_of_bad_interval [bad_interval] = stop_frame; 
            }

            number_of_samples_in_bad_interval [bad_interval] =  stop_sample_of_bad_interval [bad_interval] - start_sample_of_bad_interval [bad_interval];
        }

        

    #define SEARCH_RANGE_IN_TRANSFORM_LENGTH    4

        search_range_in_samples= SEARCH_RANGE_IN_TRANSFORM_LENGTH * Nf;

        for (bad_interval= 0; bad_interval< number_of_bad_intervals; bad_interval++) {
            float  *ref = (float *) safe_malloc ( (2 * search_range_in_samples + number_of_samples_in_bad_interval [bad_interval]) * sizeof (float));
            float  *deg = (float *) safe_malloc ( (2 * search_range_in_samples + number_of_samples_in_bad_interval [bad_interval]) * sizeof (float));
            int        i;
            float    best_correlation;
            int        delay_in_samples;

            for (i = 0; i < search_range_in_samples; i++) {
                ref[i] = 0.0f;
            }
            for (i = 0; i < number_of_samples_in_bad_interval [bad_interval]; i++) {
                ref [search_range_in_samples + i] = ref_info-> data [start_sample_of_bad_interval [bad_interval] + i];
            }
            for (i = 0; i < search_range_in_samples; i++) {
                ref [search_range_in_samples + number_of_samples_in_bad_interval [bad_interval] + i] = 0.0f;
            }
        
            for (i = 0; 
                 i < 2 * search_range_in_samples + number_of_samples_in_bad_interval [bad_interval];
                 i++) {
                
                int j = start_sample_of_bad_interval [bad_interval] - search_range_in_samples + i;
                int nn = maxNsamples - SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000);

                if (j < SEARCHBUFFER * Downsample) {
                    j = SEARCHBUFFER * Downsample;
                }
                if (j >= nn) {
                    j = nn - 1;
                }
                deg [i] = tweaked_deg [j]; 
            }

            delay_in_samples= compute_delay (0, 
                                             2 * search_range_in_samples + number_of_samples_in_bad_interval [bad_interval], 
                                             search_range_in_samples,
                                             ref, 
                                             deg,
                                             &best_correlation);

            delay_in_samples_in_bad_interval [bad_interval] =  delay_in_samples;

            if (best_correlation < 0.5) {
                delay_in_samples_in_bad_interval  [bad_interval] = 0;
            } 

            safe_free (ref);
            safe_free (deg);
        }

        if (number_of_bad_intervals > 0) {
            doubly_tweaked_deg = (float *) safe_malloc ((maxNsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof (float));

            for (i = 0; i < maxNsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
                doubly_tweaked_deg [i] = tweaked_deg [i];
            }
        
            for (bad_interval= 0; bad_interval< number_of_bad_intervals; bad_interval++) {
                int delay = delay_in_samples_in_bad_interval  [bad_interval];
                int i;
    
                for (i = start_sample_of_bad_interval [bad_interval]; i < stop_sample_of_bad_interval [bad_interval]; i++) {
                    float h;
                    int j = i + delay;
                    if (j < 0) {
                        j = 0;
                    }
                    if (j >= maxNsamples) {
                        j = maxNsamples - 1;
    
                    }
                    doubly_tweaked_deg [i] = h = tweaked_deg [j];        
                }
            }
    
            untweaked_deg = deg_info-> data;
            deg_info-> data = doubly_tweaked_deg;
    
            for (bad_interval= 0; bad_interval < number_of_bad_intervals; bad_interval++) {
    
                 for (frame = start_frame_of_bad_interval [bad_interval]; 
                     frame < stop_frame_of_bad_interval [bad_interval]; 
                     frame++) {
    
                     int start_sample_ref = SEARCHBUFFER * Downsample + frame * Nf / 2;
                    int start_sample_deg = start_sample_ref;
                    
                    short_term_fft (Nf, deg_info, Whanning, start_sample_deg, hz_spectrum_deg, fft_tmp);            
    
                    freq_warping (Nf / 2, hz_spectrum_deg, Nb, pitch_pow_dens_deg, frame);
                }    
    
                oldScale = 1;
                for (frame = start_frame_of_bad_interval [bad_interval]; 
                     frame < stop_frame_of_bad_interval [bad_interval]; 
                     frame++) {
                    int band;
    
                    total_audible_pow_ref = total_audible (frame, pitch_pow_dens_ref, 1);
                    total_audible_pow_deg = total_audible (frame, pitch_pow_dens_deg, 1);        
    
                    scale = (total_audible_pow_ref + (float) 5E3) / (total_audible_pow_deg + (float) 5E3);
                    
                    if (frame > 0) {
                        scale = (float) 0.2 * oldScale + (float) 0.8*scale;
                    }
                    oldScale = scale;
    
                    if (scale > (float) MAX_SCALE) scale = (float) MAX_SCALE;
    
                    if (scale < (float) MIN_SCALE) {
                        scale = (float) MIN_SCALE;            
                    }
    
                    for (band = 0; band < Nb; band++) {
                        pitch_pow_dens_deg [frame * Nb + band] *= scale;
                    }
    
                    intensity_warping_of (loudness_dens_ref, frame, pitch_pow_dens_ref); 
                    intensity_warping_of (loudness_dens_deg, frame, pitch_pow_dens_deg); 
        
                    for (band = 0; band < Nb; band++) {
                        disturbance_dens [band] = loudness_dens_deg [band] - loudness_dens_ref [band];
                    }
    
                    for (band = 0; band < Nb; band++) {
                        deadzone [band] = min (loudness_dens_deg [band], loudness_dens_ref [band]);    
                        deadzone [band] *= 0.25;
                    }
                    
                    for (band = 0; band < Nb; band++) {
                        float d = disturbance_dens [band];
                        float m = deadzone [band];
                    
                        if (d > m) {
                            disturbance_dens [band] -= m;
                        } else {
                            if (d < -m) {
                                disturbance_dens [band] += m;
                            } else {
                                disturbance_dens [band] = 0;
                            }
                        }
                    }    
    
                    frame_disturbance [frame] = min (frame_disturbance [frame] , pseudo_Lp (Nb, disturbance_dens, D_POW_F));    
    
                    multiply_with_asymmetry_factor (disturbance_dens, frame, pitch_pow_dens_ref, pitch_pow_dens_deg);
        
                    frame_disturbance_asym_add [frame] = min (frame_disturbance_asym_add [frame], pseudo_Lp (Nb, disturbance_dens, A_POW_F));    
                }
            }    
            safe_free (doubly_tweaked_deg);
            deg_info->data = untweaked_deg;
        }
    }
    

    for (frame = 0; frame <= stop_frame; frame++) {
        float h = 1;
        
        if (stop_frame + 1 > 1000) {
            long n = (maxNsamples - 2 * SEARCHBUFFER * Downsample) / (Nf / 2) - 1;
            double timeWeightFactor = (n - (float) 1000) / (float) 5500;
            if (timeWeightFactor > (float) 0.5) timeWeightFactor = (float) 0.5;
            h = (float) (((float) 1.0 - timeWeightFactor) + timeWeightFactor * (float) frame / (float) n);
        }

        time_weight [frame] = h;
    }

    for (frame = 0; frame <= stop_frame; frame++) {

        float h = (float) pow ((total_power_ref [frame] + 1E5) / 1E7, 0.04); 

        frame_disturbance [frame] /= h;
        frame_disturbance_asym_add [frame] /= h;

        if (frame_disturbance [frame] > 45) {
            frame_disturbance [frame] = 45;
        }
        if (frame_disturbance_asym_add [frame] > 45) {
            frame_disturbance_asym_add [frame] = 45;
        }            
    }
        
    d_indicator = Lpq_weight (start_frame, stop_frame, D_POW_S, D_POW_T, frame_disturbance, time_weight);    
    a_indicator = Lpq_weight (start_frame, stop_frame, A_POW_S, A_POW_T, frame_disturbance_asym_add, time_weight);       
    
    err_info-> pesq_mos = (float) (4.5 - D_WEIGHT * d_indicator - A_WEIGHT * a_indicator); 

    FFTFree();
    safe_free (fft_tmp);
    safe_free (hz_spectrum_ref);
    safe_free (hz_spectrum_deg);
    safe_free (silent);
    safe_free (pitch_pow_dens_ref);
    safe_free (pitch_pow_dens_deg);
    safe_free (frame_was_skipped);
    safe_free (avg_pitch_pow_dens_ref);
    safe_free (avg_pitch_pow_dens_deg);
    safe_free (loudness_dens_ref);
    safe_free (loudness_dens_deg);
    safe_free (deadzone);
    safe_free (disturbance_dens);
    safe_free (disturbance_dens_asym_add);
    safe_free (total_power_ref);

    safe_free (frame_is_bad);
    safe_free (smeared_frame_is_bad);
    
    safe_free (time_weight);
    safe_free (frame_disturbance);
    safe_free (frame_disturbance_asym_add);
    safe_free (tweaked_deg);

    return;
}

/* END OF FILE */

#define ITU_RESULTS_FILE          "_pesq_itu_results.txt"
#define SIMPLE_RESULTS_FILE       "_pesq_results.txt"

#ifndef TWOPI
  #define TWOPI   6.283185307179586f
#endif

unsigned long   FFTSwapInitialised = 0;
unsigned long   FFTLog2N;
unsigned long * FFTButter;
unsigned long * FFTBitSwap;
float         * FFTPhi;

long total_malloced = 0;

void *safe_malloc (unsigned long size) {
    void *result;
    total_malloced += size;
    result = malloc (size);
    if (result == NULL) {
        printf ("malloc failed!\n");
    }
    return result;
}

void safe_free (void *p) {
    free (p);
}


void DC_block( float * data, long Nsamples )
{
    float *p;
    long count;
    float facc = 0.0f;

    long ofs = SEARCHBUFFER * Downsample;

    p = data + ofs;
    for( count = (Nsamples - 2 * ofs); count > 0L; count-- )
        facc += *(p++);
    facc /= Nsamples;

    p = data + ofs;
    for( count = (Nsamples - 2 * ofs); count > 0L; count-- )
        *(p++) -= facc;

    p = data + ofs;
    for( count = 0L; count < Downsample; count++ )
       *(p++) *= (0.5f + count) / Downsample;

    p = data + Nsamples - ofs - 1L;
    for( count = 0L; count < Downsample; count++ )
       *(p--) *= (0.5f + count) / Downsample;
}

long InIIR_Nsos;
float *InIIR_Hsos;

void apply_filters( float * data, long Nsamples )
{
    IIRFilt( InIIR_Hsos, InIIR_Nsos, NULL,
             data, Nsamples + DATAPADDING_MSECS  * (Fs / 1000), NULL );
}

float interpolate (float    freq, 
                   double   filter_curve_db [][2],
                   int      number_of_points) {
    double  result;
    int     i;
    double  freqLow, freqHigh;
    double  curveLow, curveHigh;
    
    if (freq <= filter_curve_db [0][0]) {
        freqLow = filter_curve_db [0][0];
        curveLow = filter_curve_db [0][1];
        freqHigh = filter_curve_db [1][0];
        curveHigh = filter_curve_db [1][1];

        result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);
    
        return (float) result;
    }

    if (freq >= filter_curve_db [number_of_points-1][0]) {
        freqLow = filter_curve_db [number_of_points-2][0];
        curveLow = filter_curve_db [number_of_points-2][1];
        freqHigh = filter_curve_db [number_of_points-1][0];
        curveHigh = filter_curve_db [number_of_points-1][1];

        result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);
    
        return (float) result;
    }
        
    i = 1;
    freqHigh = filter_curve_db [i][0];
    while (freqHigh < freq) {
        i++;
        freqHigh = filter_curve_db [i][0];    
    }
    curveHigh = filter_curve_db [i][1];

    freqLow = filter_curve_db [i-1][0];
    curveLow = filter_curve_db [i-1][1];

    result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);

    return (float) result;
}       


void apply_filter ( float * data, long maxNsamples, int number_of_points, double filter_curve_db [][2] )
{ 
    long    n           = maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000);
    long    pow_of_2    = nextpow2 (n);
    float    *x            = (float *) safe_malloc ((pow_of_2 + 2) * sizeof (float));

    float    factorDb, factor;
    
    float   overallGainFilter = interpolate ((float) 1000, filter_curve_db, number_of_points); 
    float   freq_resolution;
    int        i;
    
    for (i = 0; i < pow_of_2 + 2; i++) {
        x [i] = 0;
    }

    for (i = 0; i < n; i++) {
        x [i] = data [i + SEARCHBUFFER * Downsample];    
    }

    RealFFT (x, pow_of_2);
    
    freq_resolution = (float) Fs / (float) pow_of_2;


    for (i = 0; i <= pow_of_2/2; i++) { 
        factorDb = interpolate (i * freq_resolution, filter_curve_db, number_of_points) - overallGainFilter;
        factor = (float) pow ((float) 10, factorDb / (float) 20); 

        x [2 * i] *= factor;       
        x [2 * i + 1] *= factor;   
    }

    RealIFFT (x, pow_of_2);

    for (i = 0; i < n; i++) {
        data [i + SEARCHBUFFER * Downsample] = x[i];    
    }

    safe_free (x);
}

void apply_VAD( SIGNAL_INFO * pinfo, float * data, float * VAD, float * logVAD )
{
    float g;
    float LevelThresh;
    float LevelNoise;
    float StDNoise;
    float LevelSig;
    float LevelMin;
    long  count;
    long  iteration;
    long  length;
    long  start;
    long  finish;
    long  Nwindows = (*pinfo).Nsamples / Downsample;

    for( count = 0L; count < Nwindows; count++ )
    {
        VAD[count] = 0.0f;
        for( iteration = 0L; iteration < Downsample; iteration++ )
        {
            g = data[count * Downsample + iteration];
            VAD[count] += (g * g);
        }
        VAD[count] /= Downsample;
    }

    LevelThresh = 0.0f;
    for( count = 0L; count < Nwindows; count++ )
        LevelThresh += VAD[count];
    LevelThresh /= Nwindows;

    LevelMin = 0.0f;
    for( count = 0L; count < Nwindows; count++ )
        if( VAD[count] > LevelMin )
            LevelMin = VAD[count];
    if( LevelMin > 0.0f )
        LevelMin *= 1.0e-4f;
    else
        LevelMin = 1.0f;
    
    for( count = 0L; count < Nwindows; count++ )
        if( VAD[count] < LevelMin )
            VAD[count] = LevelMin;

    for( iteration = 0L; iteration < 12L; iteration++ )
    {
        LevelNoise = 0.0f;
        StDNoise = 0.0f;
        length = 0L;
        for( count = 0L; count < Nwindows; count++ )
            if( VAD[count] <= LevelThresh )
            {
                LevelNoise += VAD[count];
                length++;
            }
        if( length > 0L )
        {
            LevelNoise /= length;
            for( count = 0L; count < Nwindows; count++ )
                if( VAD[count] <= LevelThresh )
                {
                    g = VAD[count] - LevelNoise;
                    StDNoise += g * g;
                }
            StDNoise = (float)sqrt(StDNoise / length);
        }

        LevelThresh = 1.001f * (LevelNoise + 2.0f * StDNoise);
    }

    LevelNoise = 0.0f;
    LevelSig = 0.0f;
    length = 0L;
    for( count = 0L; count < Nwindows; count++ )
    {
        if( VAD[count] > LevelThresh )
        {
            LevelSig += VAD[count];
            length++;
        }
        else
            LevelNoise += VAD[count];
    }
    if( length > 0L )
        LevelSig /= length;
    else
        LevelThresh = -1.0f;
    if( length < Nwindows )
        LevelNoise /= (Nwindows - length);
    else
        LevelNoise = 1.0f;

    for( count = 0L; count < Nwindows; count++ )
        if( VAD[count] <= LevelThresh )
            VAD[count] = -VAD[count];

    VAD[0] = -LevelMin;
    VAD[Nwindows-1] = -LevelMin;

    start = 0L;
    finish = 0L;
    for( count = 1; count < Nwindows; count++ )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
            start = count;
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
        {
            finish = count;
            if( (finish - start) <= MINSPEECHLGTH )
                for( iteration = start; iteration < finish; iteration++ )
                    VAD[iteration] = -VAD[iteration];
        }
    }

    if( LevelSig >= (LevelNoise * 1000.0f) )
    {
        for( count = 1; count < Nwindows; count++ )
        {
            if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
                start = count;
            if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
            {
                finish = count;
                g = 0.0f;
                for( iteration = start; iteration < finish; iteration++ )
                    g += VAD[iteration];
                if( g < 3.0f * LevelThresh * (finish - start) )
                    for( iteration = start; iteration < finish; iteration++ )
                        VAD[iteration] = -VAD[iteration];
            }
        }
    }

    start = 0L;
    finish = 0L;
    for( count = 1; count < Nwindows; count++ )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
        {
            start = count;
            if( (finish > 0L) && ((start - finish) <= JOINSPEECHLGTH) )
                for( iteration = finish; iteration < start; iteration++ )
                    VAD[iteration] = LevelMin;
        }
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
            finish = count;
    }

    start = 0L;
    for( count = 1; count < Nwindows; count++ )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
            start = count;
    }
    if( start == 0L )
    {
        for( count = 0L; count < Nwindows; count++ )
            VAD[count] = (float)fabs(VAD[count]);
        VAD[0] = -LevelMin;
        VAD[Nwindows-1] = -LevelMin;
    }

    count = 3;
    while( count < (Nwindows-2) )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-2] <= 0.0f) )
        {
            VAD[count-2] = VAD[count] * 0.1f;
            VAD[count-1] = VAD[count] * 0.3f;
            count++;
        }
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
        {
            VAD[count] = VAD[count-1] * 0.3f;
            VAD[count+1] = VAD[count-1] * 0.1f;
            count += 3;
        }
        count++;
    }

    for( count = 0L; count < Nwindows; count++ )
        if( VAD[count] < 0.0f ) VAD[count] = 0.0f;
    
    if( LevelThresh <= 0.0f )
        LevelThresh = LevelMin;
    for( count = 0L; count < Nwindows; count++ )
    {
        if( VAD[count] <= LevelThresh )
            logVAD[count] = 0.0f;
        else
            logVAD[count] = (float)log( VAD[count]/LevelThresh );
    }
}

void crude_align(
    SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info, ERROR_INFO * err_info,
    long Utt_id, float * ftmp)
{
    long  nr;
    long  nd;
    long  startr;
    long  startd;
    long  count;
    long  I_max;
    float max;
    float * ref_VAD = (*ref_info).logVAD;
    float * deg_VAD = (*deg_info).logVAD;
    float * Y;

    if( Utt_id == WHOLE_SIGNAL )
    {
        nr = (*ref_info).Nsamples / Downsample;
        nd = (*deg_info).Nsamples / Downsample;
        startr = 0L;
        startd = 0L;
    }
    else if( Utt_id == MAXNUTTERANCES )
    {
        startr = (*err_info).UttSearch_Start[MAXNUTTERANCES-1];
        startd = startr + (*err_info).Utt_DelayEst[MAXNUTTERANCES-1] / Downsample;

        if ( startd < 0L )
        {
            startr = -(*err_info).Utt_DelayEst[MAXNUTTERANCES-1] / Downsample;
            startd = 0L;
        }

        nr = (*err_info).UttSearch_End[MAXNUTTERANCES-1] - startr;
        nd = nr;

        if( startd + nd > (*deg_info).Nsamples / Downsample )
            nd = (*deg_info).Nsamples / Downsample - startd;
    }
    else
    {
        startr = (*err_info).UttSearch_Start[Utt_id];
        startd = startr + (*err_info).Crude_DelayEst / Downsample;

        if ( startd < 0L )
        {
            startr = -(*err_info).Crude_DelayEst / Downsample;
            startd = 0L;
        }

        nr = (*err_info).UttSearch_End[Utt_id] - startr;
        nd = nr;

        if( startd + nd > (*deg_info).Nsamples / Downsample )
            nd = (*deg_info).Nsamples / Downsample - startd;
    }

    Y  = ftmp;
        
    if( (nr > 1L) && (nd > 1L) )
        FFTNXCorr( ref_VAD + startr, nr, deg_VAD + startd, nd, Y );

    max = 0.0f;
    I_max = nr - 1;
    if( (nr > 1L) && (nd > 1L) )
        for( count = 0L; count < (nr+nd-1); count++ )
            if( Y[count] > max )
            {
                max = Y[count];
                I_max = count;
            }

    if( Utt_id == WHOLE_SIGNAL )
    {
        (*err_info).Crude_DelayEst = (I_max - nr + 1) * Downsample;
        (*err_info).Crude_DelayConf = 0.0f;
    }
    else if( Utt_id == MAXNUTTERANCES )
    {
        (*err_info).Utt_Delay[MAXNUTTERANCES-1] =
            (I_max - nr + 1) * Downsample + (*err_info).Utt_DelayEst[MAXNUTTERANCES-1];
    }
    else
    {
        (*err_info).Utt_DelayEst[Utt_id] =
            (I_max - nr + 1) * Downsample + (*err_info).Crude_DelayEst;
    }

    FFTFree();
}

void time_align(
    SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info, ERROR_INFO * err_info,
    long Utt_id, float * ftmp )
{
    long  count;
    long  I_max;
    float v_max;
    long  estdelay;
    long  startr;
    long  startd;
    float * X1;
    float * X2;
    float * H;
    float * Window;
    float r1, i1;
    long  kernel;
    float Hsum;

    estdelay = (*err_info).Utt_DelayEst[Utt_id];

    X1 = ftmp;
    X2 = ftmp + Align_Nfft + 2;
    H  = (ftmp + 4 + 2 * Align_Nfft);
    for( count = 0L; count < Align_Nfft; count++ )
        H[count] = 0.0f;
    Window = ftmp + 5 * Align_Nfft;

    for( count = 0L; count < Align_Nfft; count++ )
         Window[count] = (float)(0.5 * (1.0 - cos((TWOPI * count) / Align_Nfft)));

    startr = (*err_info).UttSearch_Start[Utt_id] * Downsample;
    startd = startr + estdelay;

    if ( startd < 0L )
    {
        startr = -estdelay;
        startd = 0L;
    }

    while( ((startd + Align_Nfft) <= (*deg_info).Nsamples) &&
           ((startr + Align_Nfft) <= ((*err_info).UttSearch_End[Utt_id] * Downsample)) )
    {
        for( count = 0L; count < Align_Nfft; count++ )
        {
            X1[count] = (*ref_info).data[count + startr] * Window[count];
            X2[count] = (*deg_info).data[count + startd] * Window[count];
            
        }
        RealFFT( X1, Align_Nfft );
        RealFFT( X2, Align_Nfft );

        for( count = 0L; count <= Align_Nfft / 2; count++ )
        {
            r1 = X1[count * 2]; i1 = -X1[1 + (count * 2)];
            X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
            X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
        }

        RealIFFT( X1, Align_Nfft );

        v_max = 0.0f;
        for( count = 0L; count < Align_Nfft; count++ )
        {
            r1 = (float) fabs(X1[count]);
            X1[count] = r1;
            if( r1 > v_max ) v_max = r1;
        }
        v_max *= 0.99f;
        for( count = 0L; count < Align_Nfft; count++ )
            if( X1[count] > v_max )
                H[count] += (float) pow( v_max, 0.125 );

        startr += (Align_Nfft / 4);
        startd += (Align_Nfft / 4);
    }

    Hsum = 0.0f;
    for( count = 0L; count < Align_Nfft; count++ )
    {
        Hsum += H[count];
        X1[count] = H[count];
        X2[count] = 0.0f;        
    }

    X2[0] = 1.0f;
    kernel = Align_Nfft / 64;
    for( count = 1; count < kernel; count++ )
    {
        X2[count] = 1.0f - ((float)count) / ((float)kernel);
        X2[(Align_Nfft - count)] = 1.0f - ((float)count) / ((float)kernel);
    }
    RealFFT( X1, Align_Nfft );
    RealFFT( X2, Align_Nfft );

    for( count = 0L; count <= Align_Nfft / 2; count++ )
    {
        r1 = X1[count * 2]; i1 = X1[1 + (count * 2)];
        X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
        X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
    }
    RealIFFT( X1, Align_Nfft );

    for( count = 0L; count < Align_Nfft; count++ )
    {
        if( Hsum > 0.0 )
            H[count] = (float) fabs(X1[count]) / Hsum;
        else
            H[count] = 0.0f;
    }

    v_max = 0.0f;
    I_max = 0L;
    for( count = 0L; count < Align_Nfft; count++ )
        if( H[count] > v_max )
        {
            v_max = H[count];
            I_max = count;
        }
    if( I_max >= (Align_Nfft/2) )
        I_max -= Align_Nfft;

    (*err_info).Utt_Delay[Utt_id] = estdelay + I_max;
    (*err_info).Utt_DelayConf[Utt_id] = v_max;

    FFTFree();
}

void split_align( SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info, float * ftmp,
    long Utt_Start, long Utt_SpeechStart, long Utt_SpeechEnd, long Utt_End,
    long Utt_DelayEst, float Utt_DelayConf,
    long * Best_ED1, long * Best_D1, float * Best_DC1,
    long * Best_ED2, long * Best_D2, float * Best_DC2,
    long * Best_BP )
{
    long count, bp, k;
    long Utt_Len = Utt_SpeechEnd - Utt_SpeechStart;
    long Utt_Test = MAXNUTTERANCES - 1;

    long N_BPs;
    long Utt_BPs[41];
    long Utt_ED1[41], Utt_ED2[41];
    long Utt_D1[41], Utt_D2[41];
    float Utt_DC1[41], Utt_DC2[41];

    long Delta, Step, Pad;

    long  estdelay;
    long  I_max;
    float v_max, n_max;
    long  startr;
    long  startd;
    float * X1;
    float * X2;
    float * H;
    float * Window;
    float r1, i1;
    long  kernel;
    float Hsum;

    *Best_DC1 = 0.0f;
    *Best_DC2 = 0.0f;

    X1 = ftmp;
    X2 = ftmp + 2 + Align_Nfft;
    H  = (ftmp + 4 + 2 * Align_Nfft);
    Window = ftmp + 6 + 3 * Align_Nfft;
    for( count = 0L; count < Align_Nfft; count++ )
         Window[count] = (float)(0.5 * (1.0 - cos((TWOPI * count) / Align_Nfft)));
    kernel = Align_Nfft / 64;

    Delta = Align_Nfft / (4 * Downsample);

    Step = (long) ((0.801 * Utt_Len + 40 * Delta - 1)/(40 * Delta));
    Step *= Delta;

    Pad = Utt_Len / 10;
    if( Pad < 75 ) Pad = 75;
    Utt_BPs[0] = Utt_SpeechStart + Pad;
    N_BPs = 0;
    do {
        N_BPs++;
        Utt_BPs[N_BPs] = Utt_BPs[N_BPs-1] + Step;
    } while( (Utt_BPs[N_BPs] <= (Utt_SpeechEnd - Pad)) && (N_BPs < 40) );

    if( N_BPs <= 0 ) return;  

    for( bp = 0; bp < N_BPs; bp++ )
    {
        (*err_info).Utt_DelayEst[Utt_Test] = Utt_DelayEst;
        (*err_info).UttSearch_Start[Utt_Test] = Utt_Start;
        (*err_info).UttSearch_End[Utt_Test] = Utt_BPs[bp];

        crude_align( ref_info, deg_info, err_info, MAXNUTTERANCES, ftmp);
        Utt_ED1[bp] = (*err_info).Utt_Delay[Utt_Test];

        (*err_info).Utt_DelayEst[Utt_Test] = Utt_DelayEst;
        (*err_info).UttSearch_Start[Utt_Test] = Utt_BPs[bp];
        (*err_info).UttSearch_End[Utt_Test] = Utt_End;

        crude_align( ref_info, deg_info, err_info, MAXNUTTERANCES, ftmp);
        Utt_ED2[bp] = (*err_info).Utt_Delay[Utt_Test];
    }

    for( bp = 0; bp < N_BPs; bp++ )
        Utt_DC1[bp] = -2.0f;
    while( 1 )
    {
        bp = 0;
        while( (bp < N_BPs) && (Utt_DC1[bp] > -2.0) )
            bp++;
        if( bp >= N_BPs )
            break;

        estdelay = Utt_ED1[bp];

        for( count = 0L; count < Align_Nfft; count++ )
            H[count] = 0.0f;
        Hsum = 0.0f;

        startr = Utt_Start * Downsample;
        startd = startr + estdelay;

        if ( startd < 0L )
        {
            startr = -estdelay;
            startd = 0L;
        }

        while( ((startd + Align_Nfft) <= (*deg_info).Nsamples) &&
               ((startr + Align_Nfft) <= (Utt_BPs[bp] * Downsample)) )
        {
            for( count = 0L; count < Align_Nfft; count++ )
            {
                X1[count] = (*ref_info).data[count + startr] * Window[count];
                X2[count] = (*deg_info).data[count + startd] * Window[count];                
            }
            RealFFT( X1, Align_Nfft );
            RealFFT( X2, Align_Nfft );

            for( count = 0L; count <= Align_Nfft / 2; count++ )
            {
                r1 = X1[count * 2]; i1 = -X1[1 + (count * 2)];
                X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
                X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
            }

            RealIFFT( X1, Align_Nfft );

            v_max = 0.0f;
            for( count = 0L; count < Align_Nfft; count++ )
            {
                r1 = (float) fabs(X1[count]);
                X1[count] = r1;
                if( r1 > v_max ) v_max = r1;
            }
            v_max *= 0.99f;
            n_max = (float) pow( v_max, 0.125 ) / kernel;

            for( count = 0L; count < Align_Nfft; count++ )
                if( X1[count] > v_max )
                {
                    Hsum += n_max * kernel;
                    for( k = 1-kernel; k < kernel; k++ )
                        H[(count + k + Align_Nfft) % Align_Nfft] +=
                            n_max * (kernel - (float) fabs(k));
                }

            startr += (Align_Nfft / 4);
            startd += (Align_Nfft / 4);
        }

        v_max = 0.0f;
        I_max = 0L;
        for( count = 0L; count < Align_Nfft; count++ )
            if( H[count] > v_max )
            {
                v_max = H[count];
                I_max = count;
            }
        if( I_max >= (Align_Nfft/2) )
            I_max -= Align_Nfft;

        Utt_D1[bp] = estdelay + I_max;
        if( Hsum > 0.0 )
            Utt_DC1[bp] = v_max / Hsum;
        else
            Utt_DC1[bp] = 0.0f;

        while( bp < (N_BPs - 1) )
        {
            bp++;
            if( (Utt_ED1[bp] == estdelay) && (Utt_DC1[bp] <= -2.0) )
            {
                while( ((startd + Align_Nfft) <= (*deg_info).Nsamples) &&
                       ((startr + Align_Nfft) <= (Utt_BPs[bp] * Downsample)) )
                {
                    for( count = 0L; count < Align_Nfft; count++ )
                    {
                        X1[count] = (*ref_info).data[count + startr] * Window[count];
                        X2[count] = (*deg_info).data[count + startd] * Window[count];                        
                    }
                    RealFFT( X1, Align_Nfft );
                    RealFFT( X2, Align_Nfft );

                    for( count = 0L; count <= Align_Nfft/2; count++ )
                    {
                        r1 = X1[count * 2]; i1 = -X1[1 + (count * 2)];
                        X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
                        X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
                    }

                    RealIFFT( X1, Align_Nfft );

                    v_max = 0.0f;
                    for( count = 0L; count < Align_Nfft; count++ )
                    {
                        r1 = (float) fabs(X1[count]);
                        X1[count] = r1;
                        if( r1 > v_max ) v_max = r1;
                    }
                    v_max *= 0.99f;
                    n_max = (float) pow( v_max, 0.125 ) / kernel;

                    for( count = 0L; count < Align_Nfft; count++ )
                        if( X1[count] > v_max )
                        {
                            Hsum += n_max * kernel;
                            for( k = 1-kernel; k < kernel; k++ )
                                H[(count + k + Align_Nfft) % Align_Nfft] +=
                                    n_max * (kernel - (float) fabs(k));
                        }

                    startr += (Align_Nfft / 4);
                    startd += (Align_Nfft / 4);
                }

                v_max = 0.0f;
                I_max = 0L;
                for( count = 0L; count < Align_Nfft; count++ )
                    if( H[count] > v_max )
                    {
                        v_max = H[count];
                        I_max = count;
                    }
                if( I_max >= (Align_Nfft/2) )
                    I_max -= Align_Nfft;

                Utt_D1[bp] = estdelay + I_max;
                if( Hsum > 0.0 )
                    Utt_DC1[bp] = v_max / Hsum;
                else
                    Utt_DC1[bp] = 0.0f;
            }
        }
    }

    for( bp = 0; bp < N_BPs; bp++ )
    {
        if( Utt_DC1[bp] > Utt_DelayConf )
            Utt_DC2[bp] = -2.0f;
        else
            Utt_DC2[bp] = 0.0f;
    }
    while( 1 )
    {
        bp = N_BPs - 1;
        while( (bp >= 0) && (Utt_DC2[bp] > -2.0) )
            bp--;
        if( bp < 0 )
            break;

        estdelay = Utt_ED2[bp];

        for( count = 0L; count < Align_Nfft; count++ )
            H[count] = 0.0f;
        Hsum = 0.0f;

        startr = Utt_End * Downsample - Align_Nfft;
        startd = startr + estdelay;

        if ( (startd + Align_Nfft) > (*deg_info).Nsamples )
        {
            startd = (*deg_info).Nsamples - Align_Nfft;
            startr = startd - estdelay;
        }

        while( (startd >= 0L) &&
               (startr >= (Utt_BPs[bp] * Downsample)) )
        {
            for( count = 0L; count < Align_Nfft; count++ )
            {
                X1[count] = (*ref_info).data[count + startr] * Window[count];
                X2[count] = (*deg_info).data[count + startd] * Window[count];                
            }
            RealFFT( X1, Align_Nfft );
            RealFFT( X2, Align_Nfft );

            for( count = 0L; count <= Align_Nfft/2; count++ )
            {
                r1 = X1[count * 2]; i1 = -X1[1 + (count * 2)];
                X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
                X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
            }

            RealIFFT( X1, Align_Nfft );

            v_max = 0.0f;
            for( count = 0L; count < Align_Nfft; count++ )
            {
                r1 = (float) fabs(X1[count]);
                X1[count] = r1;
                if( r1 > v_max ) v_max = r1;
            }
            v_max *= 0.99f;
            n_max = (float) pow( v_max, 0.125 ) / kernel;

            for( count = 0L; count < Align_Nfft; count++ )
                if( X1[count] > v_max )
                {
                    Hsum += n_max * kernel;
                    for( k = 1-kernel; k < kernel; k++ )
                        H[(count + k + Align_Nfft) % Align_Nfft] +=
                            n_max * (kernel - (float) fabs(k));
                }

            startr -= (Align_Nfft / 4);
            startd -= (Align_Nfft / 4);
        }

        v_max = 0.0f;
        I_max = 0L;
        for( count = 0L; count < Align_Nfft; count++ )
            if( H[count] > v_max )
            {
                v_max = H[count];
                I_max = count;
            }
        if( I_max >= (Align_Nfft/2) )
            I_max -= Align_Nfft;

        Utt_D2[bp] = estdelay + I_max;
        if( Hsum > 0.0 )
            Utt_DC2[bp] = v_max / Hsum;
        else
            Utt_DC2[bp] = 0.0f;

        while( bp > 0 )
        {
            bp--;
            if( (Utt_ED2[bp] == estdelay) && (Utt_DC2[bp] <= -2.0) )
            {
                while( (startd >= 0L) &&
                       (startr >= (Utt_BPs[bp] * Downsample)) )
                {
                    for( count = 0L; count < Align_Nfft; count++ )
                    {
                        X1[count] = (*ref_info).data[count + startr] * Window[count];
                        X2[count] = (*deg_info).data[count + startd] * Window[count];                        
                    }
                    RealFFT( X1, Align_Nfft );
                    RealFFT( X2, Align_Nfft );

                    for( count = 0L; count <= Align_Nfft / 2; count++ )
                    {
                        r1 = X1[count * 2]; i1 = -X1[1 + (count * 2)];
                        X1[count * 2] = (r1 * X2[count * 2] - i1 * X2[1 + (count * 2)]);
                        X1[1 + (count * 2)] = (r1 * X2[1 + (count * 2)] + i1 * X2[count * 2]);
                    }

                    RealIFFT( X1, Align_Nfft );

                    v_max = 0.0f;
                    for( count = 0L; count < Align_Nfft; count++ )
                    {
                        r1 = (float) fabs(X1[count]);
                        X1[count] = r1;
                        if( r1 > v_max ) v_max = r1;
                    }
                    v_max *= 0.99f;
                    n_max = (float) pow( v_max, 0.125 ) / kernel;

                    for( count = 0L; count < Align_Nfft; count++ )
                        if( X1[count] > v_max )
                        {
                            Hsum += n_max * kernel;
                            for( k = 1-kernel; k < kernel; k++ )
                                H[(count + k + Align_Nfft) % Align_Nfft] +=
                                    n_max * (kernel - (float) fabs(k));
                        }

                    startr -= (Align_Nfft / 4);
                    startd -= (Align_Nfft / 4);
                }

                v_max = 0.0f;
                I_max = 0L;
                for( count = 0L; count < Align_Nfft; count++ )
                    if( H[count] > v_max )
                    {
                        v_max = H[count];
                        I_max = count;
                    }
                if( I_max >= (Align_Nfft/2) )
                    I_max -= Align_Nfft;

                Utt_D2[bp] = estdelay + I_max;
                if( Hsum > 0.0 )
                    Utt_DC2[bp] = v_max / Hsum;
                else
                    Utt_DC2[bp] = 0.0f;
            }
        }
    }

    for( bp = 0; bp < N_BPs; bp++ )
    {
        if( (abs(Utt_D2[bp] - Utt_D1[bp]) >= Downsample) &&
            ((Utt_DC1[bp] + Utt_DC2[bp]) > ((*Best_DC1) + (*Best_DC2))) &&
            (Utt_DC1[bp] > Utt_DelayConf) && (Utt_DC2[bp] > Utt_DelayConf) )
            {
                *Best_ED1 = Utt_ED1[bp]; *Best_D1 = Utt_D1[bp]; *Best_DC1 = Utt_DC1[bp];
                *Best_ED2 = Utt_ED2[bp]; *Best_D2 = Utt_D2[bp]; *Best_DC2 = Utt_DC2[bp];
                *Best_BP = Utt_BPs[bp];
            }
    }

    FFTFree();
}

/* END OF FILE */

unsigned long nextpow2(unsigned long X)
{
    unsigned long C = 1;
    while( (C < ULONG_MAX) && (C < X) )
        C <<= 1;

    return C;
}

int ispow2(unsigned long X)
{
    unsigned long C = 1;
    while( (C < ULONG_MAX) && (C < X) )
        C <<= 1;
        
    return (C == X);
}

int intlog2(unsigned long X)
{
    return (int)floor( log( 1.0 * X ) / log( 2.0 ) + 0.5 );
}

void FFTInit(unsigned long N)
{
    unsigned long   C, L, K;
    float           Theta;
    float         * PFFTPhi;
    
    if( (FFTSwapInitialised != N) && (FFTSwapInitialised != 0) )
        FFTFree();

    if( FFTSwapInitialised == N )
    {
        return;
    }
    else
    {
        C = N;
        for( FFTLog2N = 0; C > 1; C >>= 1 )
            FFTLog2N++;

        C = 1;
        C <<= FFTLog2N;
        if( N == C )
            FFTSwapInitialised = N;

        FFTButter = (unsigned long *) safe_malloc( sizeof(unsigned long) * (N >> 1) );
        FFTBitSwap = (unsigned long *) safe_malloc( sizeof(unsigned long) * N );
        FFTPhi = (float *) safe_malloc( 2 * sizeof(float) * (N >> 1) );
    
        PFFTPhi = FFTPhi;
        for( C = 0; C < (N >> 1); C++ )
        {
            Theta = (TWOPI * C) / N;
            (*(PFFTPhi++)) = (float) cos( Theta );
            (*(PFFTPhi++)) = (float) sin( Theta );
        }
    
        FFTButter[0] = 0;
        L = 1;
        K = N >> 2;
        while( K >= 1 )
        {
            for( C = 0; C < L; C++ )
                FFTButter[C+L] = FFTButter[C] + K;
            L <<= 1;
            K >>= 1;
        }
    }
}

void FFTFree(void)
{
    if( FFTSwapInitialised != 0 )
    {
        safe_free( FFTButter );
        safe_free( FFTBitSwap );
        safe_free( FFTPhi );
        FFTSwapInitialised = 0;
    }
}

void FFT(float * x, unsigned long N)
{
    unsigned long   Cycle, C, S, NC;
    unsigned long   Step    = N >> 1;
    unsigned long   K1, K2;
    register float  R1, I1, R2, I2;
    float           ReFFTPhi, ImFFTPhi;

    if( N > 1 )
    {
        FFTInit( N );
    
        for( Cycle = 1; Cycle < N; Cycle <<= 1, Step >>= 1 )
        {
            K1 = 0;
            K2 = Step << 1;
    
            for( C = 0; C < Cycle; C++ )
            {
                NC = FFTButter[C] << 1;
                ReFFTPhi = FFTPhi[NC];
                ImFFTPhi = FFTPhi[NC+1];
                for( S = 0; S < Step; S++ )
                {
                    R1 = x[K1];
                    I1 = x[K1+1];
                    R2 = x[K2];
                    I2 = x[K2+1];
                    
                    x[K1++] = R1 + ReFFTPhi * R2 + ImFFTPhi * I2;
                    x[K1++] = I1 - ImFFTPhi * R2 + ReFFTPhi * I2;
                    x[K2++] = R1 - ReFFTPhi * R2 - ImFFTPhi * I2;
                    x[K2++] = I1 + ImFFTPhi * R2 - ReFFTPhi * I2;
                }
                K1 = K2;
                K2 = K1 + (Step << 1);
            }
        }
    
        NC = N >> 1;
        for( C = 0; C < NC; C++ )
        {
            FFTBitSwap[C] = FFTButter[C] << 1;
            FFTBitSwap[C+NC] = 1 + FFTBitSwap[C];
        }
        for( C = 0; C < N; C++ )
            if( (S = FFTBitSwap[C]) != C )
            {
                FFTBitSwap[S] = S;
                K1 = C << 1;
                K2 = S << 1;
                R1 = x[K1];
                x[K1++] = x[K2];
                x[K2++] = R1;
                R1 = x[K1];
                x[K1] = x[K2];
                x[K2] = R1;
            }
    }
}

void IFFT(float * x, unsigned long N)
{
    unsigned long   Cycle, C, S, NC;
    unsigned long   Step    = N >> 1;
    unsigned long   K1, K2;
    register float  R1, I1, R2, I2;
    float           ReFFTPhi, ImFFTPhi;

    if( N > 1 )
    {
        FFTInit( N );
    
        for( Cycle = 1; Cycle < N; Cycle <<= 1, Step >>= 1 )
        {
            K1 = 0;
            K2 = Step << 1;
    
            for( C = 0; C < Cycle; C++ )
            {
                NC = FFTButter[C] << 1;
                ReFFTPhi = FFTPhi[NC];
                ImFFTPhi = FFTPhi[NC+1];
                for( S = 0; S < Step; S++ )
                {
                    R1 = x[K1];
                    I1 = x[K1+1];
                    R2 = x[K2];
                    I2 = x[K2+1];
                    
                    x[K1++] = R1 + ReFFTPhi * R2 - ImFFTPhi * I2;
                    x[K1++] = I1 + ImFFTPhi * R2 + ReFFTPhi * I2;
                    x[K2++] = R1 - ReFFTPhi * R2 + ImFFTPhi * I2;
                    x[K2++] = I1 - ImFFTPhi * R2 - ReFFTPhi * I2;
                }
                K1 = K2;
                K2 = K1 + (Step << 1);
            }
        }
    
        NC = N >> 1;
        for( C = 0; C < NC; C++ )
        {
            FFTBitSwap[C] = FFTButter[C] << 1;
            FFTBitSwap[C+NC] = 1 + FFTBitSwap[C];
        }
        for( C = 0; C < N; C++ )
            if( (S = FFTBitSwap[C]) != C )
            {
                FFTBitSwap[S] = S;
                K1 = C << 1;
                K2 = S << 1;
                R1 = x[K1];
                x[K1++] = x[K2];
                x[K2++] = R1;
                R1 = x[K1];
                x[K1] = x[K2];
                x[K2] = R1;
            }
    
        NC = N << 1;
        for( C = 0; C < NC; )
            x[C++] /= N;
    }
}

void RealFFT(float *x, unsigned long N) 
{
    float            *y;
    unsigned long    i;

    y = (float *) safe_malloc (2 * N * sizeof (float));

    for (i = 0; i < N; i++) {
        y [2 * i] = x [i];
        y [2 * i + 1] = 0.0f;
    }

    FFT (y, N);

    for (i = 0; i <= N / 2; i++) {
        x [2 * i] = y [2 * i];
        x [2 * i + 1] = y [2 * i + 1];
    }    
    
    safe_free (y);
}

void RealIFFT(float *x, unsigned long N)
{

    float            *y;
    unsigned long    i;

    y = (float *) safe_malloc (2 * N * sizeof (float));

    for (i = 0; i <= N / 2; i++) {
        y [2 * i] = x [2 * i];
        y [2 * i + 1] = x [2 * i + 1];
    }    
    for (i = N / 2 + 1; i < N; i++) {
        int j = N - i;
        y [2 * i] = x [2 * j];
        y [2 * i + 1] = -x [2 * j + 1];
    }    
    
    IFFT (y, N);
    
    for (i = 0; i < N; i++) {
        x [i] = y [2 * i];
    }

    safe_free (y);
}

unsigned long FFTNXCorr(
  float * x1, unsigned long n1,
  float * x2, unsigned long n2,
  float * y )
{
    register float  r1, i1;
    float         * tmp1;
    float         * tmp2;
    long            C, D, Nx, Ny;

    Nx = nextpow2( max(n1, n2) );
    tmp1 = (float *) safe_malloc(sizeof(float) * (2 * Nx + 2));
    tmp2 = (float *) safe_malloc(sizeof(float) * (2 * Nx + 2));

    for( C = n1 - 1; C >= 0; C-- )
    {
        tmp1[C] = *(x1++);
    }
    for( C = n1; C < 2 * Nx; C++ )
        tmp1[C] = 0.0;

    RealFFT( tmp1, 2*Nx );
    
    for( C = 0; C < (long) n2; C++ )
    {
        tmp2[C] = x2[C];
    }
    for( C = n2; C < 2 * Nx; C++ )
        tmp2[C] = 0.0;
    
    RealFFT( tmp2, 2*Nx );

    for( C = 0; C <= Nx; C++ )
    {
        D = C << 1; r1 = tmp1[D]; i1 = tmp1[1 + D];
        tmp1[D] = r1 * tmp2[D] - i1 * tmp2[1 + D];
        tmp1[1 + D] = r1 * tmp2[1 + D] + i1 * tmp2[D];
    }

    RealIFFT( tmp1, 2*Nx );
    Ny = n1 + n2 - 1;
    for( C = 0; C < Ny; C++ )
        y[C] = tmp1[C];
    
    safe_free( tmp1 );
    safe_free( tmp2 );

    return Ny;
}

void IIRsos(
    float * x, unsigned long Nx,
    float b0, float b1, float b2, float a1, float a2,
    float * tz1, float * tz2 )
{
    register float z0;
    register float z1;
    register float z2;

    if( tz1 == NULL ) z1 = 0.0f; else z1 = *tz1;
    if( tz2 == NULL ) z2 = 0.0f; else z2 = *tz2;
    
    if( (a1 != 0.0f) || (a2 != 0.0f) )
    {
        if( (b1 != 0.0f) || (b2 != 0.0f) )
        {
            while( (Nx) > 0 )
            {
                Nx--;
                                z0 = (*x) - a1 * z1 - a2 * z2;
                *(x++) = b0 * z0 + b1 * z1 + b2 * z2;
                z2 = z1;
                z1 = z0;
            }
        }
        else
        {
            if( b0 != 1.0f )
            {
                while( (Nx) > 0 )
                                {
                                        Nx--;
                    z0 = (*x) - a1 * z1 - a2 * z2;
                    *(x++) = b0 * z0;
                    z2 = z1;
                    z1 = z0;
                }
            }
            else
            {
                while( (Nx) > 0 )
                {
                    Nx--;
                    z0 = (*x) - a1 * z1 - a2 * z2;
                    *(x++) = z0;
                    z2 = z1;
                    z1 = z0;
                }
            }
        }
    }
    else
    {
        if( (b1 != 0.0f) || (b2 != 0.0f) )
        {
            while( (Nx) > 0 )
            {
                Nx--;
                z0 = (*x);
                *(x++) = b0 * z0 + b1 * z1 + b2 * z2;
                z2 = z1;
                z1 = z0;
            }
        }
        else
        {
            if( b0 != 1.0f )
            {
                while( (Nx) > 0 )
                {
                    Nx--;
                    *x = b0 * (*x);
                    x++;
                }
            }
        }
    }

    if( tz1 != NULL ) (*tz1) = z1;
    if( tz2 != NULL ) (*tz2) = z2;
}

void IIRFilt(
    float * h, unsigned long Nsos, float * z,
    float * x, unsigned long Nx, float * y )
{
    unsigned long C;

    if( y == NULL )
        y = x;
    else
    {
        for( C = 0; C < Nx; C++ )
            y[C] = x[C];
    }


    for( C = 0; C < Nsos; C++ )
    {
        if( z != NULL )
        {
            IIRsos( y, Nx, h[0], h[1], h[2], h[3], h[4], z, z+1 );
            z += 2;
        }
        else
            IIRsos( y, Nx, h[0], h[1], h[2], h[3], h[4], NULL, NULL );
        h += 5;
    }
}
/* END OF FILE */

double maiin (int argc, const char *argv []);
void usage (void);
void pesq_measure (SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info, long * Error_Flag, char ** Error_Type);

void usage (void) {
    printf ("Usage:\n");
    printf (" PESQ HELP               Displays this text\n");
    printf (" PESQ [options] ref deg [smos] [cond]\n");
    printf (" Run model on reference ref and degraded deg\n");
    printf ("\n");
    printf ("Options: +8000 +16000 +swap\n");
    printf (" Sample rate - No default. Must select either +8000 or +16000.\n");
    printf (" Swap byte order - machine native format by default. Select +swap for byteswap.\n");
    printf ("\n");
    printf (" [smos] is an optional number copied to %s\n", ITU_RESULTS_FILE);
    printf (" [cond] is an optional condition number copied to %s\n", ITU_RESULTS_FILE);
    printf (" smos must always precede cond. However, both may be omitted.");
    printf ("\n");
    printf ("File names, smos, cond may not begin with a + character.\n");
    printf ("\n");
    printf ("Files with names ending .wav or .WAV are assumed to have a 44-byte header, which");
    printf (" is automatically skipped.  All other file types are assumed to have no header.\n");
}

double maiin (int argc, const char *argv []) {
    int  arg;
    int  names = 0;
    long sample_rate = -1;
    
    SIGNAL_INFO ref_info;
    SIGNAL_INFO deg_info;
    ERROR_INFO err_info;

    long Error_Flag = 0;
    char * Error_Type = "Unknown error type.";

    if (Error_Flag == 0) {
        printf("\n");

        if (argc < 3){
            usage ();
            return -1;                                                                  
        } else {

            strcpy (ref_info.path_name, "");
            ref_info.apply_swap = 0;
            strcpy (deg_info.path_name, "");
            deg_info.apply_swap = 0;
            err_info. subj_mos = 0;
            err_info. cond_nr = 0;

            for (arg = 1; arg < argc; arg++) {
                if (argv [arg] [0] == '+') {
                    if (strcmp (argv [arg], "+swap") == 0) {
                        ref_info.apply_swap = 1;
                        deg_info.apply_swap = 1;
                    } else {
                        if (strcmp (argv [arg], "+16000") == 0) {
                            sample_rate = 16000L;
                        } else {
                            if (strcmp (argv [arg], "+8000") == 0) {
                                sample_rate = 8000L;
                            } else {
                                usage ();
                                fprintf (stderr, "Invalid parameter '%s'.\n", argv [arg]);
                                return 1;
                            }
                        }
                    }
                } else {
                    switch (names) {
                        case 0: 
                            strcpy (ref_info.path_name, argv [arg]); 
                            break;
                        case 1: 
                            strcpy (deg_info.path_name, argv [arg]); 
                            break;
                        case 2: 
                            sscanf (argv [arg], "%f", &(err_info. subj_mos)); 
                            break;
                        case 3: 
                            sscanf (argv [arg], "%d", &(err_info. cond_nr)); 
                            break;
                        default:
                            usage ();
                            fprintf (stderr, "Invalid parameter '%s'.\n", argv [arg]);
                            return 1;
                    }
                    names++;
                }
            }
            sample_rate = 8000L;
            // if (sample_rate == -1) {
            //     printf ("PESQ Error. Must specify either +8000 or +16000 sample frequency option!\n");
            //     exit (1);
            // }
            
            strcpy (ref_info. file_name, ref_info. path_name);
            if (strrchr (ref_info. file_name, '\\') != NULL) {
                strcpy (ref_info. file_name, 1 + strrchr (ref_info. file_name, '\\'));
            }
            if (strrchr (ref_info. file_name, '/') != NULL) {
                strcpy (ref_info. file_name, 1 + strrchr (ref_info. file_name, '/'));
            }                

            strcpy (deg_info. file_name, deg_info. path_name);
            if (strrchr (deg_info. file_name, '\\') != NULL) {
                strcpy (deg_info. file_name, 1 + strrchr (deg_info. file_name, '\\'));
            }
            if (strrchr (deg_info. file_name, '/') != NULL) {
                strcpy (deg_info. file_name, 1 + strrchr (deg_info. file_name, '/'));
            }                

            select_rate (sample_rate, &Error_Flag, &Error_Type);
            pesq_measure (&ref_info, &deg_info, &err_info, &Error_Flag, &Error_Type);
        }
    }

    if (Error_Flag == 0) {
        printf ("\nPrediction : PESQ_MOS = %.3f\n", (double) err_info.pesq_mos);
        return (double) err_info.pesq_mos;
    } else {
        //printf ("An error of type %d ", Error_Flag);
        if (Error_Type != NULL) {
            printf (" (%s) occurred during processing.\n", Error_Type);
        } else {
            printf ("occurred during processing.\n");
        }

        return -1;
    }
}

double align_filter_dB [26] [2] = {{0.,-500},
                                 {50., -500},
                                 {100., -500},
                                 {125., -500},
                                 {160., -500},
                                 {200., -500},
                                 {250., -500},
                                 {300., -500},
                                 {350.,  0},
                                 {400.,  0},
                                 {500.,  0},
                                 {600.,  0},
                                 {630.,  0},
                                 {800.,  0},
                                 {1000., 0},
                                 {1250., 0},
                                 {1600., 0},
                                 {2000., 0},
                                 {2500., 0},
                                 {3000., 0},
                                 {3250., 0},
                                 {3500., -500},
                                 {4000., -500},
                                 {5000., -500},
                                 {6300., -500},
                                 {8000., -500}}; 


double standard_IRS_filter_dB [26] [2] = {{  0., -200},
                                         { 50., -40}, 
                                         {100., -20},
                                         {125., -12},
                                         {160.,  -6},
                                         {200.,   0},
                                         {250.,   4},
                                         {300.,   6},
                                         {350.,   8},
                                         {400.,  10},
                                         {500.,  11},
                                         {600.,  12},
                                         {700.,  12},
                                         {800.,  12},
                                         {1000., 12},
                                         {1300., 12},
                                         {1600., 12},
                                         {2000., 12},
                                         {2500., 12},
                                         {3000., 12},
                                         {3250., 12},
                                         {3500., 4},
                                         {4000., -200},
                                         {5000., -200},
                                         {6300., -200},
                                         {8000., -200}}; 


#define TARGET_AVG_POWER    1E7

void fix_power_level (SIGNAL_INFO *info, char *name, long maxNsamples) 
{
    long   n = info-> Nsamples;
    long   i;
    float *align_filtered = (float *) safe_malloc ((n + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof (float));    
    float  global_scale;
    float  power_above_300Hz;

    for (i = 0; i < n + DATAPADDING_MSECS  * (Fs / 1000); i++) {
        align_filtered [i] = info-> data [i];
    }
    apply_filter (align_filtered, info-> Nsamples, 26, align_filter_dB);

    power_above_300Hz = (float) pow_of (align_filtered, 
                                        SEARCHBUFFER * Downsample, 
                                        n - SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000),
                                        maxNsamples - 2 * SEARCHBUFFER * Downsample + DATAPADDING_MSECS  * (Fs / 1000));

    global_scale = (float) sqrt (TARGET_AVG_POWER / power_above_300Hz); 

    for (i = 0; i < n; i++) {
        info-> data [i] *= global_scale;    
    }

    safe_free (align_filtered);
}

       
void pesq_measure (SIGNAL_INFO * ref_info, SIGNAL_INFO * deg_info,
    ERROR_INFO * err_info, long * Error_Flag, char ** Error_Type)
{
    float * ftmp = NULL;

    ref_info-> data = NULL;
    ref_info-> VAD = NULL;
    ref_info-> logVAD = NULL;
    
    deg_info-> data = NULL;
    deg_info-> VAD = NULL;
    deg_info-> logVAD = NULL;
        
    if ((*Error_Flag) == 0)
    {
        printf ("Reading reference file %s...", ref_info-> path_name);

       load_src (Error_Flag, Error_Type, ref_info);
       if ((*Error_Flag) == 0)
           printf ("done.\n");
    }
    if ((*Error_Flag) == 0)
    {
        printf ("Reading degraded file %s...", deg_info-> path_name);

       load_src (Error_Flag, Error_Type, deg_info);
       if ((*Error_Flag) == 0)
           printf ("done.\n");
    }

    if (((ref_info-> Nsamples - 2 * SEARCHBUFFER * Downsample < Fs / 4) ||
         (deg_info-> Nsamples - 2 * SEARCHBUFFER * Downsample < Fs / 4)) &&
        ((*Error_Flag) == 0))
    {
        (*Error_Flag) = 2;
        (*Error_Type) = "Reference or Degraded below 1/4 second - processing stopped ";
    }

    if ((*Error_Flag) == 0)
    {
        alloc_other (ref_info, deg_info, Error_Flag, Error_Type, &ftmp);
    }

    if ((*Error_Flag) == 0)
    {   
        int     maxNsamples = max (ref_info-> Nsamples, deg_info-> Nsamples);
        float * model_ref; 
        float * model_deg; 
        long    i;
        FILE *resultsFile;

        printf (" Level normalization...\n");            
        fix_power_level (ref_info, "reference", maxNsamples);
        fix_power_level (deg_info, "degraded", maxNsamples);

        printf (" IRS filtering...\n"); 
        apply_filter (ref_info-> data, ref_info-> Nsamples, 26, standard_IRS_filter_dB);
        apply_filter (deg_info-> data, deg_info-> Nsamples, 26, standard_IRS_filter_dB);

        model_ref = (float *) safe_malloc ((ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof (float));
        model_deg = (float *) safe_malloc ((deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof (float));

        for (i = 0; i < ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
            model_ref [i] = ref_info-> data [i];
        }
    
        for (i = 0; i < deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
            model_deg [i] = deg_info-> data [i];
        }
    
        input_filter( ref_info, deg_info, ftmp );

        printf (" Variable delay compensation...\n");            
        calc_VAD (ref_info);
        calc_VAD (deg_info);
        
        crude_align (ref_info, deg_info, err_info, WHOLE_SIGNAL, ftmp);

        utterance_locate (ref_info, deg_info, err_info, ftmp);
    
        for (i = 0; i < ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
            ref_info-> data [i] = model_ref [i];
        }
    
        for (i = 0; i < deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
            deg_info-> data [i] = model_deg [i];
        }

        safe_free (model_ref);
        safe_free (model_deg); 
    
        if ((*Error_Flag) == 0) {
            if (ref_info-> Nsamples < deg_info-> Nsamples) {
                float *new_ref = (float *) safe_malloc((deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof(float));
                long  i;
                for (i = 0; i < ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
                    new_ref [i] = ref_info-> data [i];
                }
                for (i = ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); 
                     i < deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
                    new_ref [i] = 0.0f;
                }
                safe_free (ref_info-> data);
                ref_info-> data = new_ref;
                new_ref = NULL;
            } else {
                if (ref_info-> Nsamples > deg_info-> Nsamples) {
                    float *new_deg = (float *) safe_malloc((ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000)) * sizeof(float));
                    long  i;
                    for (i = 0; i < deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
                        new_deg [i] = deg_info-> data [i];
                    }
                    for (i = deg_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); 
                         i < ref_info-> Nsamples + DATAPADDING_MSECS  * (Fs / 1000); i++) {
                        new_deg [i] = 0.0f;
                    }
                    safe_free (deg_info-> data);
                    deg_info-> data = new_deg;
                    new_deg = NULL;
                }
            }
        }        

        printf (" Acoustic model processing...\n");    
        pesq_psychoacoustic_model (ref_info, deg_info, err_info, ftmp);
    
        safe_free (ref_info-> data);
        safe_free (ref_info-> VAD);
        safe_free (ref_info-> logVAD);
        safe_free (deg_info-> data);
        safe_free (deg_info-> VAD);
        safe_free (deg_info-> logVAD);
        safe_free (ftmp);

        resultsFile = fopen (ITU_RESULTS_FILE, "at");

        if (resultsFile != NULL) {
            long start, end;

            if (0 != fseek (resultsFile, 0, SEEK_SET)) {
                printf ("Could not move to start of results file %s!\n", ITU_RESULTS_FILE);
                exit (1);
            }
            start = ftell (resultsFile);

            if (0 != fseek (resultsFile, 0, SEEK_END)) {
                printf ("Could not move to end of results file %s!\n", ITU_RESULTS_FILE);
                exit (1);
            }
            end = ftell (resultsFile);

            if (start == end) {
                fprintf (resultsFile, "REFERENCE\t DEGRADED\t PESQMOS\t PESQMOS\t SUBJMOS\t COND\t SAMPLE_FREQ\t CRUDE_DELAY\n");
                fflush (resultsFile);
            }

            fprintf (resultsFile, "%s\t ", ref_info-> path_name);
            fprintf (resultsFile, "%s\t ", deg_info-> path_name);
            fprintf (resultsFile, "SQValue=%.3f\t ", err_info->pesq_mos);
            fprintf (resultsFile, "%.3f\t ", err_info->pesq_mos);
            fprintf (resultsFile, "%.3f\t ", err_info->subj_mos);
            fprintf (resultsFile, "%d\t ", err_info->cond_nr);
            //fprintf (resultsFile, "%d\t", Fs);
            fprintf (resultsFile, "%.4f\n ", (float) err_info-> Crude_DelayEst / (float) Fs); 
            
            fclose (resultsFile);
        }

        resultsFile = fopen (SIMPLE_RESULTS_FILE, "at");

        if (resultsFile != NULL) {
            long start, end;

            if (0 != fseek (resultsFile, 0, SEEK_SET)) {
                printf ("Could not move to start of results file %s!\n", SIMPLE_RESULTS_FILE);
                exit (1);
            }
            start = ftell (resultsFile);

            if (0 != fseek (resultsFile, 0, SEEK_END)) {
                printf ("Could not move to end of results file %s!\n", SIMPLE_RESULTS_FILE);
                exit (1);
            }
            end = ftell (resultsFile);

            if (start == end) {
                fprintf (resultsFile, "DEGRADED\t PESQMOS\t SUBJMOS\t COND\t SAMPLE_FREQ\t CRUDE_DELAY\n");
                fflush (resultsFile);
            }

            fprintf (resultsFile, "%s\t ", deg_info-> file_name);
            fprintf (resultsFile, "%.3f\t ", err_info->pesq_mos);
            fprintf (resultsFile, "%.3f\t ", err_info->subj_mos);
            fprintf (resultsFile, "%d\t ", err_info->cond_nr);
            //fprintf (resultsFile, "%d\t", Fs);
            fprintf (resultsFile, "%.4f\n ", (float) err_info-> Crude_DelayEst / (float) Fs); 
            
            fclose (resultsFile);
        }
    }

    return;
}

/* END OF FILE */
