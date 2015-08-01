//! \file jpl_eph.c
//! Read JPL ephemeris files.
//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2012 Adam Hincks.
//
// This file is part of AEPHEM.
//
// AEPHEM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// AEPHEM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with AEPHEM.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================
// Some code in this file is directly based on that in the jpl_eph library of 
// Project Pluto, which is distributed under the GNU General Public License.
// Jpl_eph is available at www.projectpluto.com/jpl_eph.htm.
//==============================================================================

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aephem.h"

#define AE_JPL_MAXLINE     128  //!< The maximum number of characters for input.
#define AE_JPL_TTL_MAXLINE  86  //!< The maximum length of TTL input.
//! A macro for swopping a block.
#define AE_SWOP_MACRO(A, B, TMP)   {TMP = A; A = B; B = TMP;}
#define AE_JPL_DATA_GROUP 1070  //!< The group number for JPL ASCII data.
#define AE_JPL_GROUP_WORD "GROUP" //!< The keyword for a new ASCII group.

void ae_jpl_init_common_start(struct ae_jpl_handle_t *jh);
void ae_jpl_init_common_end(struct ae_jpl_handle_t *jh);
int ae_jpl_ascii_find_group(FILE *fp, int group_num);
int ae_jpl_close_wrap(struct ae_jpl_handle_t *jh, int ret_code);
int ae_jpl_get_block(struct ae_jpl_handle_t *jh, int block_num);
int ae_jpl_read_ascii_block(struct ae_jpl_handle_t *jh, int block_num);
int ae_jpl_read_bin_block(struct ae_jpl_handle_t *jh, int block_num);
int ae_jpl_cheby(struct ae_jpl_handle_t *jh, int obj_num, double block_time, 
                 int n_coord, double *pv);
int ae_jpl_read_ascii_float(FILE *fp, double *val);
void ae_swop_32_bit_val(void *ptr);
void ae_swop_64_bit_val(void *ptr, int n);


//------------------------------------------------------------------------------
//! Initialise to use a JPL ephemeris file.
//! This function opens the file, determines whether it is binary or ASCII, and
//! parses the header information.  This information is stored in the handle \p
//! jh and is passed to subsequent JPL function calls.
//!
//! If the ephemeris is in ASCII format and has a separate header file, use
//! ae_jpl_init_ascii() instead.
//!
//! If you know the format of the ephemeris file, you can use ae_jpl_init_bin()
//! or ae_jpl_init_ascii() directly.
//!
//! \param path The location of the ephemeris file.
//! \param jh For returning an initialised handle.
//!
//! \return 0 on success; otherwise, a negative code from #ae_retcode_t; this
//!         assuming that it is a binary file, which is the last format
//!         attempted.
//------------------------------------------------------------------------------

int ae_jpl_init(const char *path, struct ae_jpl_handle_t *jh) {
  int i;

  jh->fp = NULL;

  if (ae_jpl_init_ascii(path, NULL, 0, jh)) {
    if ((i = ae_jpl_init_bin(path, jh))) {
      jh->fp = NULL;
      return i;
    }
    jh->is_bin = 1;
  }
  else
    jh->is_bin = 0;

  return 0;
}


//------------------------------------------------------------------------------
//! Initialisation common to binary and ASCII at the beginning of the function.
//!
//! \param jh The JPL handle.
//------------------------------------------------------------------------------

void ae_jpl_init_common_start(struct ae_jpl_handle_t *jh) {
  jh->n_const = 0;
  jh->max_cheby = 0;
  jh->const_name = NULL;
  jh->const_val = NULL;
  jh->fp = NULL;
  jh->block = NULL;
  jh->pc = NULL;
  jh->vc = NULL;
}


//------------------------------------------------------------------------------
//! Initialisation common to binary and ASCII at the end of the function.
//!
//! \param jh The JPL handle.
//------------------------------------------------------------------------------

void ae_jpl_init_common_end(struct ae_jpl_handle_t *jh) {
  int i;

  // Set some values.
  jh->last_block_num = -1;
  jh->sun_pv_t = -1.0;
  jh->last_x = -1000.0;
  jh->ascii_len_block = -1;

  // Set the maximum number of Chebyshev coefficients possible for this file,
  // to initialize position and velocity Chebyshev coefficient arrays during
  // Chebyshev interpolation.
  jh->max_cheby = 0;
  for (i = 0; i < 13; i++) {
    if (jh->ipt[i][1] > jh->max_cheby)
      jh->max_cheby = jh->ipt[i][1];
  }

  // Allocate caches.
  jh->block = (double *)malloc(jh->rec_size * sizeof(double));
  jh->pc = (double *)malloc(jh->max_cheby * sizeof(double));
  jh->vc = (double *)malloc(jh->max_cheby * sizeof(double));

  return;
}

//------------------------------------------------------------------------------
//! Read a JPL ephemeris ASCII header.
//!
//! \param path The location of the ephemeris file.
//! \param header_path Set to NULL, unless the header file is separate (as it
//!                    might be an ASCII ephermeris).
//! \param search_dates Since JPL often splits its ASCII ephermerides into
//!                     multiple files, the date range specified in the header
//!                     might not correspond to the file being used.  To search
//!                     through the file to get the dates, set this to 1
//!                     (recommended). To trust the header, set to 0.
//! \param jh For returning an initialised handle, on success.
//!
//! \return 0 on success; on error, a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_jpl_init_ascii(const char *path, const char *header_path,
                      int search_dates, struct ae_jpl_handle_t *jh) {
  char curr_line[AE_JPL_MAXLINE + 1], *h_path;
  int i, j, k, this_block, this_n_coeff;
  double tmp[3];

  ae_jpl_init_common_start(jh);

  if (header_path == NULL)
    h_path = (char *)path;
  else
    h_path = (char *)header_path;
  
  if ((jh->fp = fopen(h_path, "r")) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_PATH);

  // First header line: KSIZE= NNNN NCOEFF= NNNN. We only need NCOEFF.
  if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (sscanf(curr_line, "%*s %d %*s %d", &k, &jh->n_coeff) != 2)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER);

  // Receive size is simply the number of coefficients times the size of a
  // double.
  jh->rec_size = jh->n_coeff * sizeof(double);

  // A consistency check.
  if (k != 2 * jh->n_coeff)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER_KSIZE);

  // GROUP 1010: Title of ephemeris (DE/LE number, start JD, end JD).
  // We merely check that it is here, without saving any information.
  if ((i = ae_jpl_ascii_find_group(jh->fp, 1010)))
    return ae_jpl_close_wrap(jh, i);

  // GROUP 1030: Start and End JD, timestep (in JD) per block.
  if ((i = ae_jpl_ascii_find_group(jh->fp, 1030)))
    return ae_jpl_close_wrap(jh, i);
  if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (sscanf(curr_line, "%lf %lf %lf", &jh->start, &jh->end, &jh->step) != 3)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER_RANGE);

  // GROUP 1040: Constant names.
  if ((i = ae_jpl_ascii_find_group(jh->fp, 1040)))
    return ae_jpl_close_wrap(jh, i);
  if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);

  // Allocate for constants.
  jh->n_const = atoi(curr_line);
  jh->const_name = (char **)malloc(jh->n_const * sizeof(char *));
  jh->const_val = (double *)malloc(jh->n_const * sizeof(double));
  for (i = 0; i < jh->n_const; i++)
    jh->const_name[i] = (char *)malloc((AE_JPL_CONST_LEN + 1) * sizeof(char));

  // Now read the constant names, 10 per line, each 6 characters long preceded
  // by 2 blanks.  Pad names with blanks to make 6 characters.
  for (i = 0; i < jh->n_const;) {
    if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
      return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);

    if ((j = strlen(curr_line)) < 81) {  // Pad end with blanks for copying.
      while (j < 81) 
        curr_line[j] = ' ';
    }
    for (j = 0; j < 10 && i < jh->n_const; j++, i++) {
      strncpy(jh->const_name[i], &curr_line[2 + j * 8], AE_JPL_CONST_LEN);
      jh->const_name[i][AE_JPL_CONST_LEN] = '\0';
    }
  }

  // GROUP 1041: Constant values.
  if ((i = ae_jpl_ascii_find_group(jh->fp, 1041)))
    return ae_jpl_close_wrap(jh, i);
  if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (atoi(curr_line) != jh->n_const)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER_NUM);

  // Now read constant values, 3 per line, 26 characters each.
  for (i = 0; i < jh->n_const; i += 3) {
    if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
      return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);

    // Switch D's to E's in exponent values.
    for (j = 0; j < strlen(curr_line); j++) {
      if (tolower(curr_line[j]) == 'd') 
        curr_line[j] = 'E';   // Exponent is 'E'.
    }

    // Scan in values. Be careful in case we are on the last line and the number
    // of constants is not a multiple of three.
    if (i + 3 < jh->n_const) {
      if (sscanf(curr_line, "%lf %lf %lf", &jh->const_val[i], 
                          &jh->const_val[i + 1], &jh->const_val[i + 2]) != 3)
        return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER_NUM);
    }
    else {
      j = jh->n_const % 3;
      if (j == 0)
        j = 3;
      if (sscanf(curr_line, "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]) != j)
        return ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER_NUM);
      for (k = 0; k < j; k++)
        jh->const_val[i + k] = tmp[k];
    }
  }
  
  // GROUP 1050: information how to parse out the Chebyshev coefficients.
  if ((i = ae_jpl_ascii_find_group(jh->fp, 1050)))
    return ae_jpl_close_wrap(jh, i);
  for (i = 0; i < 3; i++) {
    // Read line of 13 6-digit integers.
    if (fgets(curr_line, AE_JPL_MAXLINE, jh->fp) == NULL)
      return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
    for (j = 0; j < 13; j++)
      jh->ipt[j][i] = atoi(&curr_line[j * 6]);
  }

  // The following is taken from the ephcom-1.0 code:
  // "If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
  // then ipt[i][0] should contain the value of the next available coefficient
  // number rather than 0, as per communication of Myles Standish to Paul Hardy
  // on preferred format of ephemeris headers.
  //
  // "If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
  // should contain the value of the next available coefficient number rather
  // than 0 as well, as per the same communication from Myles Standish.
  //
  // "First set j to maximum index into ipt[] that has coefficients."
  for (i = 1, j = 0; i < 12; i++) {
    if (jh->ipt[i][1] > 0 && jh->ipt[i][0] > j)
      j = i;
    // Now set j to next available index count.
    if (jh->ipt[12][1] > 0 && jh->ipt[12][0] > j)
      j = jh->ipt[12][1] + jh->ipt[12][1] * jh->ipt[12][2] * 3;
    else
      j = jh->ipt[j][0] + jh->ipt[j][1] * jh->ipt[j][2] * (j == 11 ? 2 : 3);
    for (i = 1; i < 12; i++)
      if (jh->ipt[i][0] == 0)
        jh->ipt[i][0] = j;
    if (jh->ipt[12][0] == 0)
      jh->ipt[12][0] = j;
  }

  // Grab some of the important constant values.
  if (!(jh->au = ae_jpl_get_const_val(jh, "AU")))
    jh->au = AE_AU;     // Default value.
  if (!(jh->em_ratio = ae_jpl_get_const_val(jh, "EMRAT")))
    jh->em_ratio = AE_E_M_RAT;    // Default value.
  jh->ephemeris_version = (int)ae_jpl_get_const_val(jh, "DENUM");

  // Data group.
  if ((i = ae_jpl_ascii_find_group(jh->fp, AE_JPL_DATA_GROUP)))
     return ae_jpl_close_wrap(jh, i);

  if (header_path == NULL) {
    // The header is in the same file as the data.
    if (feof(jh->fp))
      return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);

    jh->ascii_start = ftell(jh->fp);
  }
  else {
    fclose(jh->fp);

    if ((jh->fp = fopen(path, "r")) == NULL)
      return ae_jpl_close_wrap(jh, -AE_RET_BAD_PATH);

    jh->ascii_start = 0;
  }

  ae_jpl_init_common_end(jh);
  jh->is_bin = 0;

  if (search_dates) {
    // OK, dig through the file for the start and end Julian dates.
    // Find the start date. It is the first number in the first block.
    if (fscanf(jh->fp, "%d %d", &this_block, &this_n_coeff) != 2)
      return -AE_RET_BAD_JPL_CORRUPT;
    if (this_block != 1 || this_n_coeff != jh->n_coeff)
      return -AE_RET_BAD_JPL_CORRUPT;
    if ((i = ae_jpl_read_ascii_float(jh->fp, &jh->start)))
      return i;

    // Now do the dumbest thing possible and read through the whole file until
    // EOF.
    for (i = 2; fgets(curr_line, AE_JPL_MAXLINE, jh->fp) != NULL;) {
      // The beginning of a block starts with two integers.
      if (sscanf(curr_line, "%d %d", &this_block, &this_n_coeff) == 2) {
        // Check for consistency.
        if (this_block != i || this_n_coeff != jh->n_coeff)
          return -AE_RET_BAD_JPL_CORRUPT;

        // The next number is the first Julian date for the block, and the next
        // is the final Julian date.
        if ((j = ae_jpl_read_ascii_float(jh->fp, &tmp[0])))
          return j;
        if ((j = ae_jpl_read_ascii_float(jh->fp, &jh->end)))
          return j;
        i++;
      }
    }
    
    // Go back to the beginning.
    if (fseek(jh->fp, jh->ascii_start, SEEK_SET))
      return -AE_RET_UNEXPECTED_EOF;
  }

  return 0;
}


//------------------------------------------------------------------------------
//! Read a JPL Ephemeris header in binary format.
//!
//! \param path The location of the ephemeris file.
//! \param jh For returning an initialised handle, on success.
//!
//! \return 0 on success; otherwise a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_jpl_init_bin(const char *path, struct ae_jpl_handle_t *jh) {
  int i, j, version, kernel_size;
  char title[84];
  
  ae_jpl_init_common_start(jh);

  if ((jh->fp = fopen(path, "rb")) == NULL)
    return ae_jpl_close_wrap(jh, -AE_RET_BAD_PATH);

  // Start reading the header.
  if (fread(title, 84, 1, jh->fp) != 1)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (fseek(jh->fp, 2652L, SEEK_SET))
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (fread(jh, AE_JPL_HEAD_SIZE, 1, jh->fp) != 1)
    return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);

  version = atoi(title + 26);

  // In the binary file, data is stored for ipt[0...11], then the ephemeris 
  // version,  then the remaining ipt[12] data.  A little switching is required 
  // to get the correct order.
  jh->ipt[12][0] = jh->ipt[12][1];
  jh->ipt[12][1] = jh->ipt[12][2];
  jh->ipt[12][2] = jh->ephemeris_version;
  jh->ephemeris_version = version;

  jh->swop = (jh->n_const < 0 || jh->n_const > 65536L);
  if (jh->swop) {
    // Byte order is wrong for current platform.
    ae_swop_64_bit_val(&jh->start, 1);
    ae_swop_64_bit_val(&jh->end, 1);
    ae_swop_64_bit_val(&jh->step, 1);
    ae_swop_32_bit_val(&jh->n_const);
    ae_swop_64_bit_val(&jh->au, 1);
    ae_swop_64_bit_val(&jh->em_ratio, 1);
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 13; i++)
        ae_swop_32_bit_val(&jh->ipt[i][j]);
    }
  }
  
  // A sanity check:  if the earth-moon ratio is outside reasonable bounds, 
  // we must be looking at a wrong or corrupted file.
  // In DE-102,  em_ratio = 81.3007;  in DE-405/406, em_ratio = 81.30056. */
  // Those are the low and high ranges.  We'll allow some slop in case the 
  // earth/moon mass ratio changes.
  if (jh->em_ratio > 81.3008 || jh->em_ratio < 81.30055)
    ae_jpl_close_wrap(jh, -AE_RET_BAD_JPL_HEADER);

  // Once upon a time,  the kernel size was determined from the DE version.
  // This was not a terrible idea, except that it meant that when the code 
  // faced a new version, it broke.  Now we use some logic to compute the
  // kernel size.
  kernel_size = 4;
  for (i = 0; i < 13; i++)
    kernel_size += jh->ipt[i][1] * jh->ipt[i][2] * ((i == 11) ? 4 : 6);
  jh->rec_size = kernel_size * 4;
  jh->n_coeff = kernel_size / 2;

  // Read in the constant values.
  jh->const_val = (double *)malloc(jh->n_const * sizeof(double));
  fseek(jh->fp, jh->rec_size, SEEK_SET);
  if (fread(jh->const_val, sizeof(double), (size_t)jh->n_const, jh->fp) != 
      jh->n_const)
    ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
  if (jh->swop)
    ae_swop_64_bit_val(jh->const_val, jh->n_const);

  // Read in the constant names.
  fseek(jh->fp, 84 * 3, SEEK_SET); // I.e., after the 3 title lines.
  jh->const_name = (char **)malloc(jh->n_const * sizeof(char *));
  for (i = 0; i < jh->n_const; i++) {
    jh->const_name[i] = (char *)malloc((AE_JPL_CONST_LEN + 1) * sizeof(char));
    if (fread(jh->const_name[i], AE_JPL_CONST_LEN, 1, jh->fp) != 1)
      return ae_jpl_close_wrap(jh, -AE_RET_UNEXPECTED_EOF);
    jh->const_name[i][6] = '\0';
  }

  ae_jpl_init_common_end(jh);
  jh->is_bin = 1;

  return 0;
}


//------------------------------------------------------------------------------
//! An internal wrapper for closing a JPL file.
//! Not designed to be called by the user.
//!
//! \param jh The handle to be closed.
//! \param ret_code The value this function should return.
//!
//! \return Always returns \p ret_code.
//------------------------------------------------------------------------------

int ae_jpl_close_wrap(struct ae_jpl_handle_t *jh, int ret_code) {
  int i;

  if (jh->fp != NULL)
    fclose(jh->fp);
  jh->fp = NULL;

  for (i = 0; i < jh->n_const; i++)
    free(jh->const_name[i]);
  free(jh->const_name);
  free(jh->const_val);
  free(jh->block);
  free(jh->pc);
  free(jh->vc);

  return ret_code;
}

//------------------------------------------------------------------------------
//! Finish with a JPL handle.
//! Closes any open file pointers and frees any allocated arrays.
//!
//! \param jh A JPL handle.
//------------------------------------------------------------------------------

void ae_jpl_close(struct ae_jpl_handle_t *jh) {
  int dummy;

  dummy = ae_jpl_close_wrap(jh, 0);

  return;
}


//------------------------------------------------------------------------------
//! Get the value of a constant defined in the JPL ephemeris header.
//!
//! \param jh An initialised JPL handle.
//! \param const_name The name of the constant, as it appears in the header.
//!
//! \return The constant value; 0 if the constant name does not exist.
//------------------------------------------------------------------------------

double ae_jpl_get_const_val(const struct ae_jpl_handle_t *jh,
                            const char *const_name) {
  int i, j;

  for (i = 0; i < jh->n_const; i++) {
    if (strlen(const_name) < AE_JPL_CONST_LEN)
      j = strlen(const_name);
    else
      j = AE_JPL_CONST_LEN;

    if (!strncmp(jh->const_name[i], const_name, j))
      return jh->const_val[i];
  }

  return 0;
}


//------------------------------------------------------------------------------
//! Read the next floating point number from an ASCII file.
//! Annoyingly, the JPL files use 'D' instead of 'E' for their scientific
//! notation. This means that we need to fiddle around a bit to change 'D' to
//! 'E'.
//!
//! This function assumes that file stream is currently at the number to be
//! read. After it returns, the file stream will have advanced past the number.
//!
//! \param fp The file stream.
//! \param val For returning the value.
//!
//! \return 0 on success; otherwise, a negative error code from #ae_retcode_t
//!         is returned.
//------------------------------------------------------------------------------

int ae_jpl_read_ascii_float(FILE *fp, double *val) {
  char val_str[65];
  int i, len;

  if (fscanf(fp, "%64s", val_str) != 1) {
    if (feof(fp))
      return -AE_RET_UNEXPECTED_EOF;
    else
      return -AE_RET_BAD_JPL_BLOCK;
  }

  val_str[64] = '\0';
  len = strlen(val_str);

  for (i = 0; i < len; i++) {
    if (tolower(val_str[i]) == 'd')
      val_str[i] = 'e';
  }

  if (sscanf(val_str, "%lf", val) != 1)
    return -AE_RET_BAD_JPL_FLOAT;

  return 0;
}

//------------------------------------------------------------------------------
//! Read a block of data coefficients from a JPL ASCII ephemeris file.
//! 
//! \param jh The handle for this file.  (The function ae_jpl_init() should be 
//!           called first.)
//! \param block_num The block number to read. Even though in ASCII files the
//!                  block numbers start at unity, \p block_num should be
//!                  counted from 0, for compatibility with
//!                  ae_jpl_read_bin_block(). This function will add one to the
//!                  passed block number.
//!
//! \return 0 on success; otherwise, a negative error code from #ae_retcode_t
//!         is returned.
//------------------------------------------------------------------------------

int ae_jpl_read_ascii_block(struct ae_jpl_handle_t *jh, int block_num) {
  char tmpstr[AE_JPL_MAXLINE + 1];
  int i, j, n_line, n_zero, this_block, this_n_coeff, found_start;
  double val;

  // In ASCII files, the block number counting starts at 1.
  block_num++;

  if (jh->ascii_len_block < 0) {
    // Figure out the length of the block. Since ought to be at the beginning of
    // the data, if this part fails, then the file is corrupt, and the function
    // returns an error.
    if (fseek(jh->fp, jh->ascii_start, SEEK_SET))
      return -AE_RET_UNEXPECTED_EOF;
    if (fscanf(jh->fp, "%d %d", &this_block, &this_n_coeff) != 2)
      return -AE_RET_BAD_JPL_BLOCK;
    if (this_block != 1 || this_n_coeff != jh->n_coeff)
      return -AE_RET_BAD_JPL_BLOCK;

    // There are three coefficients per line. The ASCII file is padded with 
    // zeros if n_coeff isn't divisible by three.
    if ((n_zero = 3 - (jh->n_coeff % 3)) == 3)
      n_zero = 0;

    // Go through the coefficients, plus any dummy zero padding at the end of
    // the block.
    for (i = 0; i < jh->n_coeff + n_zero; i++) {
      if ((j = ae_jpl_read_ascii_float(jh->fp, &val)))
        return j;
    }
    
    // Make sure we are on a new line.
    do {
      j = fgetc(jh->fp);
    } while (j == 10 || j == 13);
    if (fseek(jh->fp, -1, SEEK_CUR))
      return -AE_RET_UNEXPECTED_EOF;

    // OK, now we should be at the beginning of the second block.
    if ((i = ftell(jh->fp)) < 0)
      return -AE_RET_UNEXPECTED_EOF;

    // Double check that we are really at the beginning of the second block.
    if (fscanf(jh->fp, "%d %d", &this_block, &this_n_coeff) != 2)
      return -AE_RET_BAD_JPL_BLOCK;
    if (this_block != 2 || this_n_coeff != jh->n_coeff)
      return -AE_RET_BAD_JPL_BLOCK;

    // Record the block length.
    jh->ascii_len_block = i - jh->ascii_start;
  }

  // First attempt:  try to seek directly to where the block ought to start.
  found_start = 0;
  i = jh->ascii_start + (block_num - 1) * jh->ascii_len_block;
  if (!fseek(jh->fp, i, SEEK_SET)) {
    this_block = 0;
    this_n_coeff = 0;
    if (fscanf(jh->fp, "%d %d", &this_block, &this_n_coeff) == 2) {
      if (this_block == block_num && this_n_coeff == jh->n_coeff)
        found_start = 1;
    }
  }

  // Second attempt. Brute force.
  if (!found_start) {
    if (fseek(jh->fp, jh->ascii_start, SEEK_SET))
      return -AE_RET_UNEXPECTED_EOF;

    for (i = 1; i != block_num; ) {
      if (fgets(tmpstr, 128, jh->fp) == NULL)
        return -AE_RET_UNEXPECTED_EOF;
      if (sscanf(tmpstr, "%d %d", &this_block, &this_n_coeff) != 2)
        return -AE_RET_BAD_JPL_CORRUPT;

      // We must insist that the block numbers increase sequentially and
      // that the number of coefficients is correct.
      if (this_block != i || this_n_coeff != jh->n_coeff)
        return -AE_RET_BAD_JPL_CORRUPT;

      // Read in the coefficient lines.
      n_line = jh->n_coeff / 3;
      if (jh->n_coeff % 3)
        n_line++;
      for (j = 0; j < n_line; j++) {
        if (fgets(tmpstr, AE_JPL_MAXLINE, jh->fp) == NULL)
          return -AE_RET_UNEXPECTED_EOF;
      }
      i++;
    }
    if (fscanf(jh->fp, "%d %d", &this_block, &this_n_coeff) != 2)
      return -AE_RET_BAD_JPL_CORRUPT;
    if (this_block != i || this_n_coeff != jh->n_coeff)
        return -AE_RET_BAD_JPL_CORRUPT;
  }

  // Finally we're all synced up. Read the values into the buffer.
  for (i = 0; i < jh->n_coeff; i++) {
    if ((j = ae_jpl_read_ascii_float(jh->fp, &jh->block[i])))
      return j;
  }
  
  return 0;
}


//------------------------------------------------------------------------------
//! Get a data block from a JPL file.
//! Reads, as appropriate, either from a binary or ASCII file and caches result.
//! If requested block_num is already cached, then no read is done.
//!
//! \param jh The handle for this file.  (The function ae_jpl_init() should be 
//!           called first.)
//! \param block_num The data block number, starting at 0.
//!
//! \return The number of coefficients in the block, or 0 on EOF.  If there is 
//!         an error, a negative error code from #ae_retcode_t is returned.
//------------------------------------------------------------------------------

int ae_jpl_get_block(struct ae_jpl_handle_t *jh, int block_num) {
  // Don't read from disc unless we have to.
  if (block_num == jh->last_block_num)
    return 0;

  // Read from disc.
  if (jh->is_bin)
    return ae_jpl_read_bin_block(jh, block_num);
  else
    return ae_jpl_read_ascii_block(jh, block_num);
}


//------------------------------------------------------------------------------
//! Read a JPL Ephemeris data block in binary format.
//! The block number ranges from 0 on up (starting at first data block, after 
//! the 2 header blocks).  
//!
//! \param jh The handle for this file.  (The function ae_jpl_init() should be 
//!           called first.)
//! \param block_num The data block number, starting at 0.
//!
//! \return On success, 0.  If there is an error, a negative error code from 
//!         #ae_retcode_t is returned.
//------------------------------------------------------------------------------

int ae_jpl_read_bin_block(struct ae_jpl_handle_t *jh, int block_num) {
  long byte;

  // Add two to block_num to account for header.
  //byte = (block_num + 2) * jh->n_coeff * jh->rec_size;
  byte = (block_num + 2) * jh->rec_size;

  // Grab the data.
  if (fseek(jh->fp, byte, SEEK_SET))
    return -AE_RET_UNEXPECTED_EOF;
  if (fread(jh->block, sizeof(double), (size_t)jh->n_coeff, jh->fp) !=
        (size_t)jh->n_coeff) {
    if (feof(jh->fp))
      return -AE_RET_UNEXPECTED_EOF;
    else
      return -AE_RET_READ_ERROR;
  }

  // Swop 'er up, if necessary.
  if (jh->swop)
    ae_swop_64_bit_val(jh->block, jh->n_coeff);

  // Store the block number.
  jh->last_block_num = block_num;

  return 0;
}


//------------------------------------------------------------------------------
//! Find a group header.
//! This function searches for the group from the beginning of the file, until
//! it either finds the requested group, or it reaches the data record group 
//! (#AE_JPL_DATA_GROUP). In the latter case, the function fails. Otherwise, it
//! returns successfully with the file synced to the beginning of the group.
//!
//! Group header lines are 12 characters long, of the form "GROUP   NNNN". We
//! allow more flexibility and only require that they be of the form "%s %d".
//!
//! \param fp The file pointer to read.
//! \param group_num The group number to search for.
//!
//! \return 0 on success; otherwise a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_jpl_ascii_find_group(FILE *fp, int group_num) {
  char curr_line[AE_JPL_MAXLINE + 1], group_word[6];
  int this_group_num;

  rewind(fp);

  // Find the group.
  this_group_num = 0;
  while (fgets(curr_line, AE_JPL_MAXLINE, fp) != NULL) {
    if (sscanf(curr_line, "%5s %d", group_word, &this_group_num) == 2) {
      if (!strcmp(group_word, AE_JPL_GROUP_WORD)) {
        if (this_group_num == group_num) {
          this_group_num = 1;
          break;
        }
        else if (this_group_num == AE_JPL_DATA_GROUP) {
          this_group_num = 0;
          break;
        }
      }
    }
  }

  if (feof(fp) || !this_group_num)
    return -AE_RET_BAD_JPL_HEADER_GROUP;

  // Skip any blank lines.
  while (fgets(curr_line, AE_JPL_MAXLINE, fp) != NULL) {
    if (strlen(curr_line))
      break;
  }

  // Allow EOF only if it is the data group (because the data could be in
  // another file).
  if (feof(fp) && group_num != AE_JPL_DATA_GROUP)
    return -AE_RET_UNEXPECTED_EOF;

  return 0;
}


//------------------------------------------------------------------------------
//! Get object positions from a JPL ephemeris file.
//!
//! \param jh A handle initialised by ae_jpl_init().
//! \param jd_tt The Julian date in TT.
//! \param obj_num The object for which to get coordinates.  If a planetary 
//!                ephemeris file is being used, one can use the predefined
//!                #ae_ss_bodies_t constants for clarity.
//! \param r For returning the rectangular position coordinates, in AU; for
//!          nutations and librations, units are seconds of arc.
//! \param v For returning the rectangular velocity, in AU per day (or seconds
//!          of arc per day for nutations or librations)  Set to NULL if you do
//!          not require this information (though setting to NULL does not make
//!          this routine faster).
//! \param is_planetary If 0, assume these are not planetary ephemerides and
//!                     barycentric coordinates are returned.  Otherwise, it 
//!                     will be assumed that a Solar System ephemeris file is 
//!                     being used and that the object numbers are as in 
//!                     #ae_ss_bodies_t; heliocentric coordinates will be
//!                     returned.
//!
//! \return 0 on success; on failure, a negative code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_jpl_get_coords(struct ae_jpl_handle_t *jh, double jd_tt, int obj_num, 
                      double r[], double v[], char is_planetary) {
  int i, n, block_num;
  double obj[6], embary[6], moon_emb[6], block_time, conversion;

  // So maybe it's _slightly wasteful not to allow the exact starting and ending
  // times, but this way we don't need to worry about making exceptions for
  // these times.
  if (jd_tt <= jh->start || jd_tt >= jh->end)
    return -AE_RET_BAD_JPL_DATE;

  // Get the data block that contains coefficients for desired date.
  block_num = (int)((jd_tt - jh->start) / jh->step);
  if ((i = ae_jpl_get_block(jh, block_num)))
    return i;
  block_time = jd_tt - jh->block[0];

  if (is_planetary) {
    // Get the sun so that we can return things in heliocentric coordinates.
    // But skip the calculation if it is already cached.
    if (jd_tt != jh->sun_pv_t) {
      if ((i = ae_jpl_cheby(jh, AE_SS_SUN, block_time, 3, jh->sun_pv)))
        return i;
      jh->sun_pv_t = jd_tt;
    }

    // Find the barycentric coordinates.  There are a couple of special cases.
    switch (obj_num) {
      case AE_SS_EARTH: case AE_SS_MOON:
        // Get the earth-moon barycentre and the position of the moon relative
        // to this barycentre.
        if ((i = ae_jpl_cheby(jh, AE_SS_EMBARY, block_time, 3, embary)))
          return i;
        if ((i = ae_jpl_cheby(jh, AE_SS_MOON_EMB, block_time, 3, moon_emb)))
          return i;
        for (i = 0; i < 6; i++) {
          switch (obj_num) {
            case AE_SS_EARTH:
              obj[i] = embary[i] - moon_emb[i] / (1.0 + jh->em_ratio);
              break;
            case AE_SS_MOON:
              obj[i] = moon_emb[i] * jh->em_ratio / (1.0 + jh->em_ratio) +
                       embary[i];
              break;
            default:
              fprintf(stderr, "Should never have reached here.  Fix the "
                              "code!\n");
              exit(0);
          }
        }
        break;

      case AE_SS_SSBARY:
        // In barycentric coordinates, the barycentre is at the origin.
        for (i = 0; i < 6; i++)
          obj[i] = 0;
        break;

      case AE_SS_NUTATION:
        // Nutation only has two coordinates.
        if ((i = ae_jpl_cheby(jh, AE_SS_NUTATION, block_time, 2, obj)))
          return i;
        break;

      default:
        if (obj_num < 0 || obj_num >= AE_N_SS_BODIES_JPL)
          return -AE_RET_BAD_OBJ_INDEX;
        if ((i = ae_jpl_cheby(jh, obj_num, block_time, 3, obj)))
          return i;
        break;
    }

    // Convert to heliocentric coordinates.
    if (obj_num != AE_SS_NUTATION && obj_num != AE_SS_LIBRATION) {
      for (i = 0; i < 6; i++)
        obj[i] -= jh->sun_pv[i];
    }
  }
  else {
    if ((i = ae_jpl_cheby(jh, obj_num, block_time, 3, obj)))
      return i;
  }

  // Seperate positions und velocities and convert to AU or arcseconds.
  if (is_planetary && obj_num == AE_SS_NUTATION) {
    conversion = AE_RTS;
    n = 2;
  }
  else if (is_planetary && obj_num == AE_SS_LIBRATION) {
    conversion = AE_RTS;
    n = 3;
  }
  else {
    conversion = 1.0 / AE_AU;
    n = 3;
  }

  for (i = 0; i < n; i++)
    r[i] = obj[i] * conversion;
  if (v != NULL) {
    for (i = 0; i < n; i++)
      v[i] = obj[i + n] * conversion;
  }

  return 0;
}


//------------------------------------------------------------------------------
//! Interpolate at a point using Chebyshev coefficients.
//!
//! \param jh A handle initialised with ae_jpl_init().
//! \param obj_num The object index.
//! \param block_time The time, in days, from the beginning of the block at
//!                   which to interpolate.
//! \param n_coord The number of coordinates for each of position and velocity.
//! \param *pv Array for returning position (first three indices) and velocity
//!            (last three indices).
//!
//! \return 0 on success; on failure, a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_jpl_cheby(struct ae_jpl_handle_t *jh, int obj_num, double block_time, 
                 int n_coord, double *pv) {
  int i, j, sub_interval, offset, n_coeff;
  double span, sub_time, x, *y, two_x;

  n_coeff = jh->ipt[obj_num][1];
  span = jh->step / jh->ipt[obj_num][2];    // Days / subinterval.
  sub_interval = (int)(block_time / span);
  offset = jh->ipt[obj_num][0] - 1 + n_coord * jh->ipt[obj_num][1] * 
           sub_interval;
  sub_time = block_time - sub_interval * jh->step / jh->ipt[obj_num][2];
    
  // Divide days in this subblock by total days in subblock to get interval 
  // [0,1].  The right part of the expression will evaluate to a whole number:
  // subinterval lengths are all integer multiples of days in a block (all 
  // powers of 2).
  x = sub_time / span;
  x = (x + x) - 1.0;
  if (x < -1.0 || x > 1.0)
    return -AE_RET_BAD_JPL_CHEBY_TIME;
  y = &(jh->block[offset]);
  
  // Only calculate Chebyshev coefficients if they're not cached.
  if (x != jh->last_x) {
    jh->last_x = x;
    two_x = x + x;   // For Chebyshev recursion.

    // Initialize position polynomial coefficients.
    jh->pc[0] = 1.0;    // Chebyshev T[0](x) = 1.
    jh->pc[1] = x;      // Chebyshev T[1](x) = x.
    for (i = 2; i < jh->max_cheby; i++) {
      jh->pc[i] = two_x * jh->pc[i - 1] - jh->pc[i - 2];
      
      // Resolve bug with gcc generating -0.0 (also makes the smallest 
      // represented number equal to zero).
      if (jh->pc[i] * jh->pc[i] == 0.0)
        jh->pc[i] = 0.0;
    }
   
    // Initialize derivative polynomial coefficients
    jh->vc[0] = 0.0;          // d(1)/dx = 0 
    jh->vc[1] = 1.0;          // d(x)/dx = 1
    jh->vc[2] = two_x + two_x;  // d(2x^2 - 1)/dx = 4x
    for (i = 3; i < jh->max_cheby; i++)
      jh->vc[i] = two_x * jh->vc[i - 1] + jh->pc[i - 1] + jh->pc[i - 1] - 
                          jh->vc[i - 2];
  }

  // Interpolate to get position and velocity for each component.
  for (i = 0; i < n_coord; i++) { 
    pv[i] = 0.0;
    pv[n_coord + i] = 0.0;
    for (j = n_coeff - 1; j >= 0; j--) {
      pv[i] += jh->pc[j] * y[i * n_coeff + j];
      pv[n_coord + i] += jh->vc[j] * y[i * n_coeff + j];
    }
    pv[n_coord + i] *= 2.0 / span;
  }

  return 0;

  return -AE_RET_UNDER_CONSTRUCTION;
}


//------------------------------------------------------------------------------
//! Byte-swop a 32-bit value in-place.
//!
//! \param ptr The 32-bit value to be swopped.
//------------------------------------------------------------------------------

void ae_swop_32_bit_val(void *ptr) {
  char *tptr, tchar;

  tptr = (char *)ptr;

  AE_SWOP_MACRO(tptr[0], tptr[3], tchar);
  AE_SWOP_MACRO(tptr[1], tptr[2], tchar);

  return;
}


//------------------------------------------------------------------------------
//! Byte-swop a list of 64-bit values in-place.
//!
//! \param ptr The 64-bit values to be swopped.
//! \param n The number of values in \p ptr.
//------------------------------------------------------------------------------

void ae_swop_64_bit_val(void *ptr, int n) {
  char *tptr, tchar;
   
  tptr = (char *)ptr;

  while (n--) {
    AE_SWOP_MACRO( tptr[0], tptr[7], tchar);
    AE_SWOP_MACRO( tptr[1], tptr[6], tchar);
    AE_SWOP_MACRO( tptr[2], tptr[5], tchar);
    AE_SWOP_MACRO( tptr[3], tptr[4], tchar);
    tptr += 8;
  }
}
