#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 3
  sp_rates[2] = -(fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 4
  sp_rates[3] = 2.0 * (fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 3
  sp_rates[2] -= (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[3] += (fwd_rates[1] - rev_rates[1]);
  //sp 6
  sp_rates[5] = -(fwd_rates[1] - rev_rates[1]);
  //sp 8
  sp_rates[7] = (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 3
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[3] += (fwd_rates[2] - rev_rates[2]);
  //sp 7
  sp_rates[6] = (fwd_rates[2] - rev_rates[2]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 5
  sp_rates[4] = (fwd_rates[3] - rev_rates[3]) * pres_mod[1];
  //sp 6
  sp_rates[5] -= 2.0 * (fwd_rates[3] - rev_rates[3]) * pres_mod[1];

  //rxn 4
  //sp 4
  sp_rates[3] -= (fwd_rates[4] - rev_rates[4]);
  //sp 5
  sp_rates[4] -= (fwd_rates[4] - rev_rates[4]);
  //sp 6
  sp_rates[5] += (fwd_rates[4] - rev_rates[4]);
  //sp 8
  sp_rates[7] += (fwd_rates[4] - rev_rates[4]);

  //rxn 5
  //sp 4
  sp_rates[3] -= (fwd_rates[5] - rev_rates[5]) * pres_mod[2];
  //sp 7
  sp_rates[6] += (fwd_rates[5] - rev_rates[5]) * pres_mod[2];
  //sp 8
  sp_rates[7] -= (fwd_rates[5] - rev_rates[5]) * pres_mod[2];

  //rxn 6
  //sp 6
  sp_rates[5] -= (fwd_rates[6] - rev_rates[6]);
  //sp 7
  sp_rates[6] -= (fwd_rates[6] - rev_rates[6]);
  //sp 8
  sp_rates[7] += 2.0 * (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 4
  sp_rates[3] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[3];
  //sp 6
  sp_rates[5] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[3];
  //sp 8
  sp_rates[7] += (fwd_rates[7] - rev_rates[7]) * pres_mod[3];

  //rxn 8
  //sp 9
  sp_rates[8] = -(fwd_rates[8] - rev_rates[8]) * pres_mod[4];
  //sp 8
  sp_rates[7] += 2.0 * (fwd_rates[8] - rev_rates[8]) * pres_mod[4];

  //rxn 9
  //sp 9
  sp_rates[8] -= (fwd_rates[9] - rev_rates[9]);
  //sp 4
  sp_rates[3] -= (fwd_rates[9] - rev_rates[9]);
  //sp 7
  sp_rates[6] += (fwd_rates[9] - rev_rates[9]);
  //sp 8
  sp_rates[7] += (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 9
  sp_rates[8] -= (fwd_rates[10] - rev_rates[10]);
  //sp 10
  sp_rates[9] = (fwd_rates[10] - rev_rates[10]);
  //sp 3
  sp_rates[2] += (fwd_rates[10] - rev_rates[10]);
  //sp 4
  sp_rates[3] -= (fwd_rates[10] - rev_rates[10]);

  //rxn 11
  //sp 9
  sp_rates[8] -= (fwd_rates[11] - rev_rates[11]);
  //sp 10
  sp_rates[9] += (fwd_rates[11] - rev_rates[11]);
  //sp 6
  sp_rates[5] -= (fwd_rates[11] - rev_rates[11]);
  //sp 8
  sp_rates[7] += (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 9
  sp_rates[8] -= (fwd_rates[12] - rev_rates[12]);
  //sp 10
  sp_rates[9] += (fwd_rates[12] - rev_rates[12]);
  //sp 7
  sp_rates[6] += (fwd_rates[12] - rev_rates[12]);
  //sp 8
  sp_rates[7] -= (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 9
  sp_rates[8] -= (fwd_rates[13] - rev_rates[13]);
  //sp 10
  sp_rates[9] += (fwd_rates[13] - rev_rates[13]);
  //sp 7
  sp_rates[6] += (fwd_rates[13] - rev_rates[13]);
  //sp 8
  sp_rates[7] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 10
  sp_rates[9] -= (fwd_rates[14] - rev_rates[14]);
  //sp 4
  sp_rates[3] -= (fwd_rates[14] - rev_rates[14]);
  //sp 8
  sp_rates[7] += 2.0 * (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 10
  sp_rates[9] -= (fwd_rates[15] - rev_rates[15]);
  //sp 3
  sp_rates[2] += (fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[4] += (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 10
  sp_rates[9] -= (fwd_rates[16] - rev_rates[16]);
  //sp 5
  sp_rates[4] += (fwd_rates[16] - rev_rates[16]);
  //sp 6
  sp_rates[5] -= (fwd_rates[16] - rev_rates[16]);
  //sp 8
  sp_rates[7] += (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 10
  sp_rates[9] -= (fwd_rates[17] - rev_rates[17]);
  //sp 5
  sp_rates[4] += (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[6] += (fwd_rates[17] - rev_rates[17]);
  //sp 8
  sp_rates[7] -= (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 10
  sp_rates[9] -= (fwd_rates[18] - rev_rates[18]);
  //sp 5
  sp_rates[4] += (fwd_rates[18] - rev_rates[18]);
  //sp 7
  sp_rates[6] += (fwd_rates[18] - rev_rates[18]);
  //sp 8
  sp_rates[7] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 9
  sp_rates[8] += (fwd_rates[19] - rev_rates[19]);
  //sp 10
  sp_rates[9] -= 2.0 * (fwd_rates[19] - rev_rates[19]);
  //sp 5
  sp_rates[4] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 9
  sp_rates[8] += (fwd_rates[20] - rev_rates[20]);
  //sp 10
  sp_rates[9] -= 2.0 * (fwd_rates[20] - rev_rates[20]);
  //sp 5
  sp_rates[4] += (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 10
  sp_rates[9] += (fwd_rates[21] - rev_rates[21]) * pres_mod[5];
  //sp 4
  sp_rates[3] -= (fwd_rates[21] - rev_rates[21]) * pres_mod[5];
  //sp 5
  sp_rates[4] -= (fwd_rates[21] - rev_rates[21]) * pres_mod[5];

  //rxn 22
  //sp 11
  sp_rates[10] = -(fwd_rates[22] - rev_rates[22]) * pres_mod[6];
  //sp 12
  sp_rates[11] = (fwd_rates[22] - rev_rates[22]) * pres_mod[6];
  //sp 6
  sp_rates[5] -= (fwd_rates[22] - rev_rates[22]) * pres_mod[6];

  //rxn 23
  //sp 12
  sp_rates[11] += (fwd_rates[23] - rev_rates[23]);
  //sp 11
  sp_rates[10] -= (fwd_rates[23] - rev_rates[23]);
  //sp 4
  sp_rates[3] += (fwd_rates[23] - rev_rates[23]);
  //sp 8
  sp_rates[7] -= (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 12
  sp_rates[11] += (fwd_rates[24] - rev_rates[24]);
  //sp 11
  sp_rates[10] -= (fwd_rates[24] - rev_rates[24]);
  //sp 4
  sp_rates[3] += (fwd_rates[24] - rev_rates[24]);
  //sp 8
  sp_rates[7] -= (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 10
  sp_rates[9] -= (fwd_rates[25] - rev_rates[25]);
  //sp 11
  sp_rates[10] -= (fwd_rates[25] - rev_rates[25]);
  //sp 12
  sp_rates[11] += (fwd_rates[25] - rev_rates[25]);
  //sp 8
  sp_rates[7] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 11
  sp_rates[10] -= (fwd_rates[26] - rev_rates[26]);
  //sp 12
  sp_rates[11] += (fwd_rates[26] - rev_rates[26]);
  //sp 5
  sp_rates[4] -= (fwd_rates[26] - rev_rates[26]);
  //sp 6
  sp_rates[5] += (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 12
  sp_rates[11] -= (fwd_rates[27] - rev_rates[27]);
  //sp 4
  sp_rates[3] -= (fwd_rates[27] - rev_rates[27]);
  //sp 29
  sp_rates[28] = (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 11
  sp_rates[10] += (fwd_rates[28] - rev_rates[28]);
  //sp 13
  sp_rates[12] = -(fwd_rates[28] - rev_rates[28]);
  //sp 8
  sp_rates[7] += (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 12
  sp_rates[11] += (fwd_rates[29] - rev_rates[29]);
  //sp 4
  sp_rates[3] += (fwd_rates[29] - rev_rates[29]);
  //sp 13
  sp_rates[12] -= (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 4
  sp_rates[3] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[7];
  //sp 14
  sp_rates[13] = (fwd_rates[30] - rev_rates[30]) * pres_mod[7];
  //sp 15
  sp_rates[14] = -(fwd_rates[30] - rev_rates[30]) * pres_mod[7];

  //rxn 31
  //sp 3
  sp_rates[2] += (fwd_rates[31] - rev_rates[31]);
  //sp 4
  sp_rates[3] -= (fwd_rates[31] - rev_rates[31]);
  //sp 14
  sp_rates[13] -= (fwd_rates[31] - rev_rates[31]);
  //sp 15
  sp_rates[14] += (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 14
  sp_rates[13] -= (fwd_rates[32] - rev_rates[32]);
  //sp 6
  sp_rates[5] -= (fwd_rates[32] - rev_rates[32]);
  //sp 15
  sp_rates[14] += (fwd_rates[32] - rev_rates[32]);
  //sp 8
  sp_rates[7] += (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 15
  sp_rates[14] += (fwd_rates[33] - rev_rates[33]);
  //sp 14
  sp_rates[13] -= (fwd_rates[33] - rev_rates[33]);
  //sp 7
  sp_rates[6] += (fwd_rates[33] - rev_rates[33]);
  //sp 8
  sp_rates[7] -= (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 9
  sp_rates[8] += (fwd_rates[34] - rev_rates[34]);
  //sp 10
  sp_rates[9] -= (fwd_rates[34] - rev_rates[34]);
  //sp 14
  sp_rates[13] -= (fwd_rates[34] - rev_rates[34]);
  //sp 15
  sp_rates[14] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 20
  sp_rates[19] = (fwd_rates[35] - rev_rates[35]);
  //sp 21
  sp_rates[20] = -(fwd_rates[35] - rev_rates[35]);
  //sp 14
  sp_rates[13] -= (fwd_rates[35] - rev_rates[35]);
  //sp 15
  sp_rates[14] += (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 10
  sp_rates[9] -= (fwd_rates[36] - rev_rates[36]);
  //sp 5
  sp_rates[4] += (fwd_rates[36] - rev_rates[36]);
  //sp 14
  sp_rates[13] += (fwd_rates[36] - rev_rates[36]);
  //sp 15
  sp_rates[14] -= (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 14
  sp_rates[13] -= (fwd_rates[37] - rev_rates[37]);
  //sp 15
  sp_rates[14] += 2.0 * (fwd_rates[37] - rev_rates[37]);
  //sp 16
  sp_rates[15] = -(fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 17
  sp_rates[16] = -(fwd_rates[38] - rev_rates[38]);
  //sp 16
  sp_rates[15] += (fwd_rates[38] - rev_rates[38]);

  //rxn 39
  //sp 17
  sp_rates[16] -= (fwd_rates[39] - rev_rates[39]);
  //sp 16
  sp_rates[15] += (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 17
  sp_rates[16] -= (fwd_rates[40] - rev_rates[40]);
  //sp 16
  sp_rates[15] += (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 17
  sp_rates[16] -= (fwd_rates[41] - rev_rates[41]);
  //sp 16
  sp_rates[15] += (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 17
  sp_rates[16] -= (fwd_rates[42] - rev_rates[42]);
  //sp 16
  sp_rates[15] += (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 4
  sp_rates[3] += fwd_rates[43];
  //sp 5
  sp_rates[4] -= fwd_rates[43];
  //sp 8
  sp_rates[7] += fwd_rates[43];
  //sp 11
  sp_rates[10] += fwd_rates[43];
  //sp 17
  sp_rates[16] -= fwd_rates[43];

  //rxn 44
  //sp 17
  sp_rates[16] -= (fwd_rates[44] - rev_rates[43]);
  //sp 11
  sp_rates[10] += (fwd_rates[44] - rev_rates[43]);
  //sp 5
  sp_rates[4] -= (fwd_rates[44] - rev_rates[43]);
  //sp 7
  sp_rates[6] += (fwd_rates[44] - rev_rates[43]);

  //rxn 45
  //sp 17
  sp_rates[16] -= (fwd_rates[45] - rev_rates[44]);
  //sp 3
  sp_rates[2] += (fwd_rates[45] - rev_rates[44]);
  //sp 11
  sp_rates[10] += (fwd_rates[45] - rev_rates[44]);
  //sp 6
  sp_rates[5] -= (fwd_rates[45] - rev_rates[44]);

  //rxn 46
  //sp 17
  sp_rates[16] -= (fwd_rates[46] - rev_rates[45]);
  //sp 26
  sp_rates[25] = (fwd_rates[46] - rev_rates[45]);
  //sp 4
  sp_rates[3] += (fwd_rates[46] - rev_rates[45]);
  //sp 6
  sp_rates[5] -= (fwd_rates[46] - rev_rates[45]);

  //rxn 47
  //sp 17
  sp_rates[16] -= (fwd_rates[47] - rev_rates[46]);
  //sp 3
  sp_rates[2] -= (fwd_rates[47] - rev_rates[46]);
  //sp 4
  sp_rates[3] += (fwd_rates[47] - rev_rates[46]);
  //sp 15
  sp_rates[14] += (fwd_rates[47] - rev_rates[46]);

  //rxn 48
  //sp 17
  sp_rates[16] -= (fwd_rates[48] - rev_rates[47]);
  //sp 19
  sp_rates[18] = (fwd_rates[48] - rev_rates[47]);
  //sp 3
  sp_rates[2] += (fwd_rates[48] - rev_rates[47]);
  //sp 4
  sp_rates[3] -= (fwd_rates[48] - rev_rates[47]);

  //rxn 49
  //sp 17
  sp_rates[16] -= (fwd_rates[49] - rev_rates[48]);
  //sp 25
  sp_rates[24] = (fwd_rates[49] - rev_rates[48]);
  //sp 4
  sp_rates[3] += (fwd_rates[49] - rev_rates[48]);
  //sp 8
  sp_rates[7] -= (fwd_rates[49] - rev_rates[48]);

  //rxn 50
  //sp 17
  sp_rates[16] -= (fwd_rates[50] - rev_rates[49]);
  //sp 25
  sp_rates[24] += (fwd_rates[50] - rev_rates[49]);
  //sp 11
  sp_rates[10] += (fwd_rates[50] - rev_rates[49]);
  //sp 12
  sp_rates[11] -= (fwd_rates[50] - rev_rates[49]);

  //rxn 51
  //sp 4
  sp_rates[3] -= (fwd_rates[51] - rev_rates[50]) * pres_mod[8];
  //sp 15
  sp_rates[14] += (fwd_rates[51] - rev_rates[50]) * pres_mod[8];
  //sp 16
  sp_rates[15] -= (fwd_rates[51] - rev_rates[50]) * pres_mod[8];

  //rxn 52
  //sp 26
  sp_rates[25] += (fwd_rates[52] - rev_rates[51]);
  //sp 5
  sp_rates[4] -= (fwd_rates[52] - rev_rates[51]);
  //sp 8
  sp_rates[7] += (fwd_rates[52] - rev_rates[51]);
  //sp 16
  sp_rates[15] -= (fwd_rates[52] - rev_rates[51]);

  //rxn 53
  //sp 12
  sp_rates[11] += fwd_rates[53];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[53];
  //sp 5
  sp_rates[4] -= fwd_rates[53];
  //sp 16
  sp_rates[15] -= fwd_rates[53];

  //rxn 54
  //sp 11
  sp_rates[10] += fwd_rates[54];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[54];
  //sp 6
  sp_rates[5] -= fwd_rates[54];
  //sp 16
  sp_rates[15] -= fwd_rates[54];

  //rxn 55
  //sp 19
  sp_rates[18] += (fwd_rates[55] - rev_rates[52]);
  //sp 3
  sp_rates[2] += (fwd_rates[55] - rev_rates[52]);
  //sp 4
  sp_rates[3] -= (fwd_rates[55] - rev_rates[52]);
  //sp 16
  sp_rates[15] -= (fwd_rates[55] - rev_rates[52]);

  //rxn 56
  //sp 19
  sp_rates[18] += (fwd_rates[56] - rev_rates[53]);
  //sp 16
  sp_rates[15] -= (fwd_rates[56] - rev_rates[53]);
  //sp 7
  sp_rates[6] += (fwd_rates[56] - rev_rates[53]);
  //sp 8
  sp_rates[7] -= (fwd_rates[56] - rev_rates[53]);

  //rxn 57
  //sp 26
  sp_rates[25] += (fwd_rates[57] - rev_rates[54]);
  //sp 19
  sp_rates[18] -= (fwd_rates[57] - rev_rates[54]);
  //sp 5
  sp_rates[4] -= (fwd_rates[57] - rev_rates[54]);
  //sp 6
  sp_rates[5] += (fwd_rates[57] - rev_rates[54]);

  //rxn 58
  //sp 11
  sp_rates[10] += (fwd_rates[58] - rev_rates[55]);
  //sp 19
  sp_rates[18] -= (fwd_rates[58] - rev_rates[55]);
  //sp 4
  sp_rates[3] += (fwd_rates[58] - rev_rates[55]);
  //sp 6
  sp_rates[5] -= (fwd_rates[58] - rev_rates[55]);

  //rxn 59
  //sp 3
  sp_rates[2] += (fwd_rates[59] - rev_rates[56]);
  //sp 18
  sp_rates[17] = (fwd_rates[59] - rev_rates[56]);
  //sp 19
  sp_rates[18] -= (fwd_rates[59] - rev_rates[56]);
  //sp 4
  sp_rates[3] -= (fwd_rates[59] - rev_rates[56]);

  //rxn 60
  //sp 26
  sp_rates[25] += (fwd_rates[60] - rev_rates[57]);
  //sp 19
  sp_rates[18] -= (fwd_rates[60] - rev_rates[57]);
  //sp 4
  sp_rates[3] += (fwd_rates[60] - rev_rates[57]);
  //sp 8
  sp_rates[7] -= (fwd_rates[60] - rev_rates[57]);

  //rxn 61
  //sp 25
  sp_rates[24] += (fwd_rates[61] - rev_rates[58]);
  //sp 19
  sp_rates[18] -= (fwd_rates[61] - rev_rates[58]);
  //sp 4
  sp_rates[3] += (fwd_rates[61] - rev_rates[58]);
  //sp 7
  sp_rates[6] -= (fwd_rates[61] - rev_rates[58]);

  //rxn 62
  //sp 11
  sp_rates[10] += (fwd_rates[62] - rev_rates[59]);
  //sp 26
  sp_rates[25] += (fwd_rates[62] - rev_rates[59]);
  //sp 19
  sp_rates[18] -= (fwd_rates[62] - rev_rates[59]);
  //sp 12
  sp_rates[11] -= (fwd_rates[62] - rev_rates[59]);

  //rxn 63
  //sp 18
  sp_rates[17] -= (fwd_rates[63] - rev_rates[60]);
  //sp 11
  sp_rates[10] += (fwd_rates[63] - rev_rates[60]);
  //sp 4
  sp_rates[3] += (fwd_rates[63] - rev_rates[60]);
  //sp 8
  sp_rates[7] -= (fwd_rates[63] - rev_rates[60]);

  //rxn 64
  //sp 18
  sp_rates[17] -= (fwd_rates[64] - rev_rates[61]);
  //sp 11
  sp_rates[10] += (fwd_rates[64] - rev_rates[61]);
  //sp 5
  sp_rates[4] -= (fwd_rates[64] - rev_rates[61]);
  //sp 6
  sp_rates[5] += (fwd_rates[64] - rev_rates[61]);

  //rxn 65
  //sp 5
  sp_rates[4] -= (fwd_rates[65] - rev_rates[62]) * pres_mod[9];
  //sp 21
  sp_rates[20] += (fwd_rates[65] - rev_rates[62]) * pres_mod[9];
  //sp 15
  sp_rates[14] -= (fwd_rates[65] - rev_rates[62]) * pres_mod[9];

  //rxn 66
  //sp 5
  sp_rates[4] -= (fwd_rates[66] - rev_rates[63]);
  //sp 6
  sp_rates[5] += (fwd_rates[66] - rev_rates[63]);
  //sp 15
  sp_rates[14] -= (fwd_rates[66] - rev_rates[63]);
  //sp 23
  sp_rates[22] = (fwd_rates[66] - rev_rates[63]);

  //rxn 67
  //sp 25
  sp_rates[24] += (fwd_rates[67] - rev_rates[64]);
  //sp 5
  sp_rates[4] -= (fwd_rates[67] - rev_rates[64]);
  //sp 15
  sp_rates[14] -= (fwd_rates[67] - rev_rates[64]);
  //sp 8
  sp_rates[7] += (fwd_rates[67] - rev_rates[64]);

  //rxn 68
  //sp 25
  sp_rates[24] += (fwd_rates[68] - rev_rates[65]);
  //sp 4
  sp_rates[3] += (fwd_rates[68] - rev_rates[65]);
  //sp 6
  sp_rates[5] -= (fwd_rates[68] - rev_rates[65]);
  //sp 15
  sp_rates[14] -= (fwd_rates[68] - rev_rates[65]);

  //rxn 69
  //sp 17
  sp_rates[16] += (fwd_rates[69] - rev_rates[66]);
  //sp 7
  sp_rates[6] += (fwd_rates[69] - rev_rates[66]);
  //sp 15
  sp_rates[14] -= (fwd_rates[69] - rev_rates[66]);
  //sp 8
  sp_rates[7] -= (fwd_rates[69] - rev_rates[66]);

  //rxn 70
  //sp 25
  sp_rates[24] += (fwd_rates[70] - rev_rates[67]);
  //sp 3
  sp_rates[2] += (fwd_rates[70] - rev_rates[67]);
  //sp 15
  sp_rates[14] -= (fwd_rates[70] - rev_rates[67]);
  //sp 8
  sp_rates[7] -= (fwd_rates[70] - rev_rates[67]);

  //rxn 71
  //sp 4
  sp_rates[3] += (fwd_rates[71] - rev_rates[68]);
  //sp 24
  sp_rates[23] = (fwd_rates[71] - rev_rates[68]);
  //sp 15
  sp_rates[14] -= (fwd_rates[71] - rev_rates[68]);
  //sp 8
  sp_rates[7] -= (fwd_rates[71] - rev_rates[68]);

  //rxn 72
  //sp 4
  sp_rates[3] += (fwd_rates[72] - rev_rates[69]);
  //sp 23
  sp_rates[22] += (fwd_rates[72] - rev_rates[69]);
  //sp 15
  sp_rates[14] -= (fwd_rates[72] - rev_rates[69]);
  //sp 8
  sp_rates[7] -= (fwd_rates[72] - rev_rates[69]);

  //rxn 73
  //sp 7
  sp_rates[6] += (fwd_rates[73] - rev_rates[70]);
  //sp 16
  sp_rates[15] += (fwd_rates[73] - rev_rates[70]);
  //sp 15
  sp_rates[14] -= (fwd_rates[73] - rev_rates[70]);
  //sp 8
  sp_rates[7] -= (fwd_rates[73] - rev_rates[70]);

  //rxn 74
  //sp 10
  sp_rates[9] -= (fwd_rates[74] - rev_rates[71]);
  //sp 23
  sp_rates[22] += (fwd_rates[74] - rev_rates[71]);
  //sp 15
  sp_rates[14] -= (fwd_rates[74] - rev_rates[71]);
  //sp 8
  sp_rates[7] += (fwd_rates[74] - rev_rates[71]);

  //rxn 75
  //sp 5
  sp_rates[4] += (fwd_rates[75] - rev_rates[72]);
  //sp 21
  sp_rates[20] -= (fwd_rates[75] - rev_rates[72]);
  //sp 6
  sp_rates[5] -= (fwd_rates[75] - rev_rates[72]);
  //sp 23
  sp_rates[22] += (fwd_rates[75] - rev_rates[72]);

  //rxn 76
  //sp 4
  sp_rates[3] -= (fwd_rates[76] - rev_rates[73]);
  //sp 21
  sp_rates[20] -= (fwd_rates[76] - rev_rates[73]);
  //sp 23
  sp_rates[22] += (fwd_rates[76] - rev_rates[73]);
  //sp 8
  sp_rates[7] += (fwd_rates[76] - rev_rates[73]);

  //rxn 77
  //sp 22
  sp_rates[21] = (fwd_rates[77] - rev_rates[74]);
  //sp 21
  sp_rates[20] -= (fwd_rates[77] - rev_rates[74]);
  //sp 5
  sp_rates[4] += (fwd_rates[77] - rev_rates[74]);
  //sp 8
  sp_rates[7] -= (fwd_rates[77] - rev_rates[74]);

  //rxn 78
  //sp 10
  sp_rates[9] -= (fwd_rates[78] - rev_rates[75]);
  //sp 20
  sp_rates[19] += (fwd_rates[78] - rev_rates[75]);
  //sp 21
  sp_rates[20] -= (fwd_rates[78] - rev_rates[75]);
  //sp 5
  sp_rates[4] += (fwd_rates[78] - rev_rates[75]);

  //rxn 79
  //sp 9
  sp_rates[8] -= (fwd_rates[79] - rev_rates[76]);
  //sp 10
  sp_rates[9] += (fwd_rates[79] - rev_rates[76]);
  //sp 20
  sp_rates[19] += (fwd_rates[79] - rev_rates[76]);
  //sp 21
  sp_rates[20] -= (fwd_rates[79] - rev_rates[76]);

  //rxn 80
  //sp 21
  sp_rates[20] -= (fwd_rates[80] - rev_rates[77]);
  //sp 15
  sp_rates[14] -= (fwd_rates[80] - rev_rates[77]);
  //sp 23
  sp_rates[22] += 2.0 * (fwd_rates[80] - rev_rates[77]);

  //rxn 81
  //sp 25
  sp_rates[24] += fwd_rates[81];
  //sp 22
  sp_rates[21] += fwd_rates[81];
  //sp 21
  sp_rates[20] -= 2.0 * fwd_rates[81];
  //sp 5
  sp_rates[4] += fwd_rates[81];

  //rxn 82
  //sp 21
  sp_rates[20] -= 2.0 * fwd_rates[82];
  //sp 5
  sp_rates[4] += fwd_rates[82];
  //sp 23
  sp_rates[22] += 2.0 * fwd_rates[82];

  //rxn 83
  //sp 20
  sp_rates[19] += (fwd_rates[83] - rev_rates[78]);
  //sp 3
  sp_rates[2] -= (fwd_rates[83] - rev_rates[78]);
  //sp 4
  sp_rates[3] += (fwd_rates[83] - rev_rates[78]);
  //sp 21
  sp_rates[20] -= (fwd_rates[83] - rev_rates[78]);

  //rxn 84
  //sp 20
  sp_rates[19] -= (fwd_rates[84] - rev_rates[79]);
  //sp 23
  sp_rates[22] += (fwd_rates[84] - rev_rates[79]);
  //sp 8
  sp_rates[7] += (fwd_rates[84] - rev_rates[79]);

  //rxn 85
  //sp 22
  sp_rates[21] -= (fwd_rates[85] - rev_rates[80]) * pres_mod[10];
  //sp 15
  sp_rates[14] += (fwd_rates[85] - rev_rates[80]) * pres_mod[10];
  //sp 8
  sp_rates[7] += (fwd_rates[85] - rev_rates[80]) * pres_mod[10];

  //rxn 86
  //sp 17
  sp_rates[16] += (fwd_rates[86] - rev_rates[81]) * pres_mod[11];
  //sp 22
  sp_rates[21] -= (fwd_rates[86] - rev_rates[81]) * pres_mod[11];
  //sp 7
  sp_rates[6] += (fwd_rates[86] - rev_rates[81]) * pres_mod[11];

  //rxn 87
  //sp 4
  sp_rates[3] += (fwd_rates[87] - rev_rates[82]) * pres_mod[12];
  //sp 22
  sp_rates[21] -= (fwd_rates[87] - rev_rates[82]) * pres_mod[12];
  //sp 24
  sp_rates[23] += (fwd_rates[87] - rev_rates[82]) * pres_mod[12];

  //rxn 88
  //sp 3
  sp_rates[2] += (fwd_rates[88] - rev_rates[83]);
  //sp 4
  sp_rates[3] -= (fwd_rates[88] - rev_rates[83]);
  //sp 22
  sp_rates[21] -= (fwd_rates[88] - rev_rates[83]);
  //sp 23
  sp_rates[22] += (fwd_rates[88] - rev_rates[83]);

  //rxn 89
  //sp 3
  sp_rates[2] += (fwd_rates[89] - rev_rates[84]);
  //sp 4
  sp_rates[3] -= (fwd_rates[89] - rev_rates[84]);
  //sp 22
  sp_rates[21] -= (fwd_rates[89] - rev_rates[84]);
  //sp 24
  sp_rates[23] += (fwd_rates[89] - rev_rates[84]);

  //rxn 90
  //sp 22
  sp_rates[21] -= (fwd_rates[90] - rev_rates[85]);
  //sp 6
  sp_rates[5] -= (fwd_rates[90] - rev_rates[85]);
  //sp 23
  sp_rates[22] += (fwd_rates[90] - rev_rates[85]);
  //sp 8
  sp_rates[7] += (fwd_rates[90] - rev_rates[85]);

  //rxn 91
  //sp 22
  sp_rates[21] -= (fwd_rates[91] - rev_rates[86]);
  //sp 24
  sp_rates[23] += (fwd_rates[91] - rev_rates[86]);
  //sp 6
  sp_rates[5] -= (fwd_rates[91] - rev_rates[86]);
  //sp 8
  sp_rates[7] += (fwd_rates[91] - rev_rates[86]);

  //rxn 92
  //sp 23
  sp_rates[22] += (fwd_rates[92] - rev_rates[87]);
  //sp 22
  sp_rates[21] -= (fwd_rates[92] - rev_rates[87]);
  //sp 7
  sp_rates[6] += (fwd_rates[92] - rev_rates[87]);
  //sp 8
  sp_rates[7] -= (fwd_rates[92] - rev_rates[87]);

  //rxn 93
  //sp 24
  sp_rates[23] += (fwd_rates[93] - rev_rates[88]);
  //sp 22
  sp_rates[21] -= (fwd_rates[93] - rev_rates[88]);
  //sp 7
  sp_rates[6] += (fwd_rates[93] - rev_rates[88]);
  //sp 8
  sp_rates[7] -= (fwd_rates[93] - rev_rates[88]);

  //rxn 94
  //sp 10
  sp_rates[9] += (fwd_rates[94] - rev_rates[89]);
  //sp 5
  sp_rates[4] -= (fwd_rates[94] - rev_rates[89]);
  //sp 22
  sp_rates[21] -= (fwd_rates[94] - rev_rates[89]);
  //sp 23
  sp_rates[22] += (fwd_rates[94] - rev_rates[89]);

  //rxn 95
  //sp 10
  sp_rates[9] += (fwd_rates[95] - rev_rates[90]);
  //sp 5
  sp_rates[4] -= (fwd_rates[95] - rev_rates[90]);
  //sp 22
  sp_rates[21] -= (fwd_rates[95] - rev_rates[90]);
  //sp 24
  sp_rates[23] += (fwd_rates[95] - rev_rates[90]);

  //rxn 96
  //sp 9
  sp_rates[8] += (fwd_rates[96] - rev_rates[91]);
  //sp 10
  sp_rates[9] -= (fwd_rates[96] - rev_rates[91]);
  //sp 22
  sp_rates[21] -= (fwd_rates[96] - rev_rates[91]);
  //sp 23
  sp_rates[22] += (fwd_rates[96] - rev_rates[91]);

  //rxn 97
  //sp 9
  sp_rates[8] += (fwd_rates[97] - rev_rates[92]);
  //sp 10
  sp_rates[9] -= (fwd_rates[97] - rev_rates[92]);
  //sp 22
  sp_rates[21] -= (fwd_rates[97] - rev_rates[92]);
  //sp 24
  sp_rates[23] += (fwd_rates[97] - rev_rates[92]);

  //rxn 98
  //sp 14
  sp_rates[13] += (fwd_rates[98] - rev_rates[93]);
  //sp 22
  sp_rates[21] -= (fwd_rates[98] - rev_rates[93]);
  //sp 15
  sp_rates[14] -= (fwd_rates[98] - rev_rates[93]);
  //sp 24
  sp_rates[23] += (fwd_rates[98] - rev_rates[93]);

  //rxn 99
  //sp 14
  sp_rates[13] += (fwd_rates[99] - rev_rates[94]);
  //sp 22
  sp_rates[21] -= (fwd_rates[99] - rev_rates[94]);
  //sp 15
  sp_rates[14] -= (fwd_rates[99] - rev_rates[94]);
  //sp 23
  sp_rates[22] += (fwd_rates[99] - rev_rates[94]);

  //rxn 100
  //sp 25
  sp_rates[24] += (fwd_rates[100] - rev_rates[95]);
  //sp 26
  sp_rates[25] -= (fwd_rates[100] - rev_rates[95]);
  //sp 22
  sp_rates[21] -= (fwd_rates[100] - rev_rates[95]);
  //sp 24
  sp_rates[23] += (fwd_rates[100] - rev_rates[95]);

  //rxn 101
  //sp 23
  sp_rates[22] -= (fwd_rates[101] - rev_rates[96]);
  //sp 24
  sp_rates[23] += (fwd_rates[101] - rev_rates[96]);

  //rxn 102
  //sp 20
  sp_rates[19] += (fwd_rates[102] - rev_rates[97]);
  //sp 21
  sp_rates[20] -= (fwd_rates[102] - rev_rates[97]);
  //sp 22
  sp_rates[21] -= (fwd_rates[102] - rev_rates[97]);
  //sp 24
  sp_rates[23] += (fwd_rates[102] - rev_rates[97]);

  //rxn 103
  //sp 25
  sp_rates[24] += (fwd_rates[103] - rev_rates[98]);
  //sp 10
  sp_rates[9] += (fwd_rates[103] - rev_rates[98]);
  //sp 5
  sp_rates[4] -= (fwd_rates[103] - rev_rates[98]);
  //sp 24
  sp_rates[23] -= (fwd_rates[103] - rev_rates[98]);

  //rxn 104
  //sp 25
  sp_rates[24] += (fwd_rates[104] - rev_rates[99]);
  //sp 10
  sp_rates[9] += (fwd_rates[104] - rev_rates[99]);
  //sp 5
  sp_rates[4] -= (fwd_rates[104] - rev_rates[99]);
  //sp 24
  sp_rates[23] -= (fwd_rates[104] - rev_rates[99]);

  //rxn 105
  //sp 25
  sp_rates[24] += (fwd_rates[105] - rev_rates[100]);
  //sp 3
  sp_rates[2] += (fwd_rates[105] - rev_rates[100]);
  //sp 4
  sp_rates[3] -= (fwd_rates[105] - rev_rates[100]);
  //sp 24
  sp_rates[23] -= (fwd_rates[105] - rev_rates[100]);

  //rxn 106
  //sp 9
  sp_rates[8] += (fwd_rates[106] - rev_rates[101]);
  //sp 10
  sp_rates[9] -= (fwd_rates[106] - rev_rates[101]);
  //sp 25
  sp_rates[24] += (fwd_rates[106] - rev_rates[101]);
  //sp 24
  sp_rates[23] -= (fwd_rates[106] - rev_rates[101]);

  //rxn 107
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[107] - rev_rates[102]);
  //sp 26
  sp_rates[25] -= (fwd_rates[107] - rev_rates[102]);
  //sp 24
  sp_rates[23] -= (fwd_rates[107] - rev_rates[102]);

  //rxn 108
  //sp 26
  sp_rates[25] -= (fwd_rates[108] - rev_rates[103]);
  //sp 11
  sp_rates[10] += (fwd_rates[108] - rev_rates[103]);
  //sp 22
  sp_rates[21] += (fwd_rates[108] - rev_rates[103]);
  //sp 24
  sp_rates[23] -= (fwd_rates[108] - rev_rates[103]);

  //rxn 109
  //sp 25
  sp_rates[24] += (fwd_rates[109] - rev_rates[104]);
  //sp 22
  sp_rates[21] += (fwd_rates[109] - rev_rates[104]);
  //sp 23
  sp_rates[22] -= (fwd_rates[109] - rev_rates[104]);
  //sp 24
  sp_rates[23] -= (fwd_rates[109] - rev_rates[104]);

  //rxn 110
  //sp 25
  sp_rates[24] += (fwd_rates[110] - rev_rates[105]);
  //sp 24
  sp_rates[23] -= (fwd_rates[110] - rev_rates[105]);
  //sp 7
  sp_rates[6] += (fwd_rates[110] - rev_rates[105]);
  //sp 8
  sp_rates[7] -= (fwd_rates[110] - rev_rates[105]);

  //rxn 111
  //sp 25
  sp_rates[24] += (fwd_rates[111] - rev_rates[106]);
  //sp 8
  sp_rates[7] += (fwd_rates[111] - rev_rates[106]);
  //sp 6
  sp_rates[5] -= (fwd_rates[111] - rev_rates[106]);
  //sp 24
  sp_rates[23] -= (fwd_rates[111] - rev_rates[106]);

  //rxn 112
  //sp 25
  sp_rates[24] += (fwd_rates[112] - rev_rates[107]);
  //sp 22
  sp_rates[21] += (fwd_rates[112] - rev_rates[107]);
  //sp 24
  sp_rates[23] -= 2.0 * (fwd_rates[112] - rev_rates[107]);

  //rxn 113
  //sp 25
  sp_rates[24] += (fwd_rates[113] - rev_rates[108]);
  //sp 10
  sp_rates[9] += (fwd_rates[113] - rev_rates[108]);
  //sp 5
  sp_rates[4] -= (fwd_rates[113] - rev_rates[108]);
  //sp 23
  sp_rates[22] -= (fwd_rates[113] - rev_rates[108]);

  //rxn 114
  //sp 25
  sp_rates[24] += (fwd_rates[114] - rev_rates[109]);
  //sp 3
  sp_rates[2] += (fwd_rates[114] - rev_rates[109]);
  //sp 4
  sp_rates[3] -= (fwd_rates[114] - rev_rates[109]);
  //sp 23
  sp_rates[22] -= (fwd_rates[114] - rev_rates[109]);

  //rxn 115
  //sp 9
  sp_rates[8] += (fwd_rates[115] - rev_rates[110]);
  //sp 10
  sp_rates[9] -= (fwd_rates[115] - rev_rates[110]);
  //sp 23
  sp_rates[22] -= (fwd_rates[115] - rev_rates[110]);
  //sp 25
  sp_rates[24] += (fwd_rates[115] - rev_rates[110]);

  //rxn 116
  //sp 25
  sp_rates[24] += (fwd_rates[116] - rev_rates[111]);
  //sp 14
  sp_rates[13] += (fwd_rates[116] - rev_rates[111]);
  //sp 15
  sp_rates[14] -= (fwd_rates[116] - rev_rates[111]);
  //sp 23
  sp_rates[22] -= (fwd_rates[116] - rev_rates[111]);

  //rxn 117
  //sp 25
  sp_rates[24] += (fwd_rates[117] - rev_rates[112]);
  //sp 22
  sp_rates[21] += (fwd_rates[117] - rev_rates[112]);
  //sp 23
  sp_rates[22] -= 2.0 * (fwd_rates[117] - rev_rates[112]);

  //rxn 118
  //sp 25
  sp_rates[24] += (fwd_rates[118] - rev_rates[113]) * pres_mod[13];
  //sp 26
  sp_rates[25] -= (fwd_rates[118] - rev_rates[113]) * pres_mod[13];
  //sp 4
  sp_rates[3] -= (fwd_rates[118] - rev_rates[113]) * pres_mod[13];

  //rxn 119
  //sp 25
  sp_rates[24] += (fwd_rates[119] - rev_rates[114]) * pres_mod[14];
  //sp 3
  sp_rates[2] -= (fwd_rates[119] - rev_rates[114]) * pres_mod[14];
  //sp 11
  sp_rates[10] -= (fwd_rates[119] - rev_rates[114]) * pres_mod[14];

  //rxn 120
  //sp 25
  sp_rates[24] -= (fwd_rates[120] - rev_rates[115]);
  //sp 10
  sp_rates[9] += (fwd_rates[120] - rev_rates[115]);
  //sp 5
  sp_rates[4] -= (fwd_rates[120] - rev_rates[115]);
  //sp 26
  sp_rates[25] += (fwd_rates[120] - rev_rates[115]);

  //rxn 121
  //sp 25
  sp_rates[24] -= (fwd_rates[121] - rev_rates[116]);
  //sp 26
  sp_rates[25] += (fwd_rates[121] - rev_rates[116]);
  //sp 6
  sp_rates[5] -= (fwd_rates[121] - rev_rates[116]);
  //sp 8
  sp_rates[7] += (fwd_rates[121] - rev_rates[116]);

  //rxn 122
  //sp 25
  sp_rates[24] -= (fwd_rates[122] - rev_rates[117]);
  //sp 26
  sp_rates[25] += (fwd_rates[122] - rev_rates[117]);
  //sp 3
  sp_rates[2] += (fwd_rates[122] - rev_rates[117]);
  //sp 4
  sp_rates[3] -= (fwd_rates[122] - rev_rates[117]);

  //rxn 123
  //sp 25
  sp_rates[24] -= (fwd_rates[123] - rev_rates[118]);
  //sp 26
  sp_rates[25] += (fwd_rates[123] - rev_rates[118]);
  //sp 7
  sp_rates[6] += (fwd_rates[123] - rev_rates[118]);
  //sp 8
  sp_rates[7] -= (fwd_rates[123] - rev_rates[118]);

  //rxn 124
  //sp 25
  sp_rates[24] -= (fwd_rates[124] - rev_rates[119]);
  //sp 10
  sp_rates[9] -= (fwd_rates[124] - rev_rates[119]);
  //sp 9
  sp_rates[8] += (fwd_rates[124] - rev_rates[119]);
  //sp 26
  sp_rates[25] += (fwd_rates[124] - rev_rates[119]);

  //rxn 125
  //sp 25
  sp_rates[24] -= (fwd_rates[125] - rev_rates[120]);
  //sp 26
  sp_rates[25] += (fwd_rates[125] - rev_rates[120]);
  //sp 14
  sp_rates[13] += (fwd_rates[125] - rev_rates[120]);
  //sp 15
  sp_rates[14] -= (fwd_rates[125] - rev_rates[120]);

  //rxn 126
  //sp 25
  sp_rates[24] -= (fwd_rates[126] - rev_rates[121]);
  //sp 26
  sp_rates[25] += (fwd_rates[126] - rev_rates[121]);
  //sp 28
  sp_rates[27] = (fwd_rates[126] - rev_rates[121]);
  //sp 29
  sp_rates[28] -= (fwd_rates[126] - rev_rates[121]);

  //rxn 127
  //sp 25
  sp_rates[24] -= (fwd_rates[127] - rev_rates[122]);
  //sp 26
  sp_rates[25] += (fwd_rates[127] - rev_rates[122]);
  //sp 22
  sp_rates[21] += (fwd_rates[127] - rev_rates[122]);
  //sp 23
  sp_rates[22] -= (fwd_rates[127] - rev_rates[122]);

  //rxn 128
  //sp 25
  sp_rates[24] -= (fwd_rates[128] - rev_rates[123]);
  //sp 26
  sp_rates[25] += (fwd_rates[128] - rev_rates[123]);
  //sp 20
  sp_rates[19] += (fwd_rates[128] - rev_rates[123]);
  //sp 21
  sp_rates[20] -= (fwd_rates[128] - rev_rates[123]);

  //rxn 129
  //sp 26
  sp_rates[25] -= (fwd_rates[129] - rev_rates[124]) * pres_mod[15];
  //sp 11
  sp_rates[10] += (fwd_rates[129] - rev_rates[124]) * pres_mod[15];
  //sp 4
  sp_rates[3] += (fwd_rates[129] - rev_rates[124]) * pres_mod[15];

  //rxn 130
  //sp 26
  sp_rates[25] -= (fwd_rates[130] - rev_rates[125]);
  //sp 11
  sp_rates[10] += (fwd_rates[130] - rev_rates[125]);
  //sp 5
  sp_rates[4] -= (fwd_rates[130] - rev_rates[125]);
  //sp 10
  sp_rates[9] += (fwd_rates[130] - rev_rates[125]);

  //rxn 131
  //sp 26
  sp_rates[25] -= (fwd_rates[131] - rev_rates[126]);
  //sp 11
  sp_rates[10] += (fwd_rates[131] - rev_rates[126]);
  //sp 6
  sp_rates[5] -= (fwd_rates[131] - rev_rates[126]);
  //sp 8
  sp_rates[7] += (fwd_rates[131] - rev_rates[126]);

  //rxn 132
  //sp 11
  sp_rates[10] += (fwd_rates[132] - rev_rates[127]);
  //sp 26
  sp_rates[25] -= (fwd_rates[132] - rev_rates[127]);
  //sp 3
  sp_rates[2] += (fwd_rates[132] - rev_rates[127]);
  //sp 4
  sp_rates[3] -= (fwd_rates[132] - rev_rates[127]);

  //rxn 133
  //sp 26
  sp_rates[25] -= (fwd_rates[133] - rev_rates[128]);
  //sp 11
  sp_rates[10] += (fwd_rates[133] - rev_rates[128]);
  //sp 7
  sp_rates[6] += (fwd_rates[133] - rev_rates[128]);
  //sp 8
  sp_rates[7] -= (fwd_rates[133] - rev_rates[128]);

  //rxn 134
  //sp 26
  sp_rates[25] -= (fwd_rates[134] - rev_rates[129]);
  //sp 11
  sp_rates[10] += (fwd_rates[134] - rev_rates[129]);
  //sp 14
  sp_rates[13] += (fwd_rates[134] - rev_rates[129]);
  //sp 15
  sp_rates[14] -= (fwd_rates[134] - rev_rates[129]);

  //rxn 135
  //sp 25
  sp_rates[24] += (fwd_rates[135] - rev_rates[130]);
  //sp 26
  sp_rates[25] -= 2.0 * (fwd_rates[135] - rev_rates[130]);
  //sp 11
  sp_rates[10] += (fwd_rates[135] - rev_rates[130]);

  //rxn 136
  //sp 12
  sp_rates[11] += (fwd_rates[136] - rev_rates[131]);
  //sp 26
  sp_rates[25] -= (fwd_rates[136] - rev_rates[131]);
  //sp 4
  sp_rates[3] += (fwd_rates[136] - rev_rates[131]);
  //sp 6
  sp_rates[5] -= (fwd_rates[136] - rev_rates[131]);

  //rxn 137
  //sp 4
  sp_rates[3] += fwd_rates[137];
  //sp 8
  sp_rates[7] += fwd_rates[137];
  //sp 10
  sp_rates[9] -= fwd_rates[137];
  //sp 12
  sp_rates[11] += fwd_rates[137];
  //sp 26
  sp_rates[25] -= fwd_rates[137];

  //rxn 138
  //sp 26
  sp_rates[25] -= 2.0 * fwd_rates[138];
  //sp 3
  sp_rates[2] += fwd_rates[138];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[138];

  //rxn 139
  //sp 25
  sp_rates[24] -= (fwd_rates[139] - rev_rates[132]) * pres_mod[16];
  //sp 4
  sp_rates[3] -= (fwd_rates[139] - rev_rates[132]) * pres_mod[16];
  //sp 24
  sp_rates[23] += (fwd_rates[139] - rev_rates[132]) * pres_mod[16];

  //rxn 140
  //sp 25
  sp_rates[24] += (fwd_rates[140] - rev_rates[133]) * pres_mod[17];
  //sp 4
  sp_rates[3] += (fwd_rates[140] - rev_rates[133]) * pres_mod[17];
  //sp 23
  sp_rates[22] -= (fwd_rates[140] - rev_rates[133]) * pres_mod[17];

  //rxn 141
  //sp 11
  sp_rates[10] += (fwd_rates[141] - rev_rates[134]);
  //sp 28
  sp_rates[27] -= (fwd_rates[141] - rev_rates[134]);
  //sp 7
  sp_rates[6] += (fwd_rates[141] - rev_rates[134]);

  //rxn 142
  //sp 12
  sp_rates[11] += (fwd_rates[142] - rev_rates[135]);
  //sp 3
  sp_rates[2] += (fwd_rates[142] - rev_rates[135]);
  //sp 28
  sp_rates[27] -= (fwd_rates[142] - rev_rates[135]);

  //rxn 143
  //sp 10
  sp_rates[9] -= (fwd_rates[143] - rev_rates[136]);
  //sp 28
  sp_rates[27] += (fwd_rates[143] - rev_rates[136]);
  //sp 29
  sp_rates[28] -= (fwd_rates[143] - rev_rates[136]);
  //sp 5
  sp_rates[4] += (fwd_rates[143] - rev_rates[136]);

  //rxn 144
  //sp 9
  sp_rates[8] -= (fwd_rates[144] - rev_rates[137]);
  //sp 10
  sp_rates[9] += (fwd_rates[144] - rev_rates[137]);
  //sp 28
  sp_rates[27] += (fwd_rates[144] - rev_rates[137]);
  //sp 29
  sp_rates[28] -= (fwd_rates[144] - rev_rates[137]);

  //rxn 145
  //sp 28
  sp_rates[27] -= fwd_rates[145];
  //sp 12
  sp_rates[11] += fwd_rates[145];
  //sp 3
  sp_rates[2] += fwd_rates[145];

  //rxn 146
  //sp 3
  sp_rates[2] += fwd_rates[146];
  //sp 4
  sp_rates[3] -= fwd_rates[146];
  //sp 8
  sp_rates[7] += fwd_rates[146];
  //sp 11
  sp_rates[10] += fwd_rates[146];
  //sp 28
  sp_rates[27] -= fwd_rates[146];

  //rxn 147
  //sp 11
  sp_rates[10] += fwd_rates[147];
  //sp 28
  sp_rates[27] -= fwd_rates[147];
  //sp 6
  sp_rates[5] -= fwd_rates[147];
  //sp 8
  sp_rates[7] += 2.0 * fwd_rates[147];

  //rxn 148
  //sp 4
  sp_rates[3] += fwd_rates[148];
  //sp 7
  sp_rates[6] += fwd_rates[148];
  //sp 8
  sp_rates[7] -= fwd_rates[148];
  //sp 12
  sp_rates[11] += fwd_rates[148];
  //sp 28
  sp_rates[27] -= fwd_rates[148];

  //rxn 149
  //sp 11
  sp_rates[10] += fwd_rates[149];
  //sp 28
  sp_rates[27] -= fwd_rates[149];
  //sp 7
  sp_rates[6] += fwd_rates[149];

  //rxn 150
  //sp 8
  sp_rates[7] += fwd_rates[150];
  //sp 11
  sp_rates[10] += fwd_rates[150];
  //sp 14
  sp_rates[13] += fwd_rates[150];
  //sp 15
  sp_rates[14] -= fwd_rates[150];
  //sp 28
  sp_rates[27] -= fwd_rates[150];

  //rxn 151
  //sp 8
  sp_rates[7] += fwd_rates[151];
  //sp 9
  sp_rates[8] += fwd_rates[151];
  //sp 10
  sp_rates[9] -= fwd_rates[151];
  //sp 11
  sp_rates[10] += fwd_rates[151];
  //sp 28
  sp_rates[27] -= fwd_rates[151];

  //rxn 152
  //sp 27
  sp_rates[26] = (fwd_rates[152] - rev_rates[138]);
  //sp 29
  sp_rates[28] -= (fwd_rates[152] - rev_rates[138]);
  //sp 8
  sp_rates[7] -= (fwd_rates[152] - rev_rates[138]);

  //rxn 153
  //sp 30
  sp_rates[29] = (fwd_rates[153] - rev_rates[139]) * pres_mod[18];
  //sp 15
  sp_rates[14] -= 2.0 * (fwd_rates[153] - rev_rates[139]) * pres_mod[18];

  //rxn 154
  //sp 4
  sp_rates[3] -= (fwd_rates[154] - rev_rates[140]) * pres_mod[19];
  //sp 30
  sp_rates[29] += (fwd_rates[154] - rev_rates[140]) * pres_mod[19];
  //sp 31
  sp_rates[30] = -(fwd_rates[154] - rev_rates[140]) * pres_mod[19];

  //rxn 155
  //sp 10
  sp_rates[9] += (fwd_rates[155] - rev_rates[141]);
  //sp 5
  sp_rates[4] -= (fwd_rates[155] - rev_rates[141]);
  //sp 30
  sp_rates[29] -= (fwd_rates[155] - rev_rates[141]);
  //sp 31
  sp_rates[30] += (fwd_rates[155] - rev_rates[141]);

  //rxn 156
  //sp 30
  sp_rates[29] -= (fwd_rates[156] - rev_rates[142]);
  //sp 6
  sp_rates[5] -= (fwd_rates[156] - rev_rates[142]);
  //sp 31
  sp_rates[30] += (fwd_rates[156] - rev_rates[142]);
  //sp 8
  sp_rates[7] += (fwd_rates[156] - rev_rates[142]);

  //rxn 157
  //sp 3
  sp_rates[2] += (fwd_rates[157] - rev_rates[143]);
  //sp 4
  sp_rates[3] -= (fwd_rates[157] - rev_rates[143]);
  //sp 30
  sp_rates[29] -= (fwd_rates[157] - rev_rates[143]);
  //sp 31
  sp_rates[30] += (fwd_rates[157] - rev_rates[143]);

  //rxn 158
  //sp 31
  sp_rates[30] += (fwd_rates[158] - rev_rates[144]);
  //sp 30
  sp_rates[29] -= (fwd_rates[158] - rev_rates[144]);
  //sp 7
  sp_rates[6] += (fwd_rates[158] - rev_rates[144]);
  //sp 8
  sp_rates[7] -= (fwd_rates[158] - rev_rates[144]);

  //rxn 159
  //sp 9
  sp_rates[8] += (fwd_rates[159] - rev_rates[145]);
  //sp 10
  sp_rates[9] -= (fwd_rates[159] - rev_rates[145]);
  //sp 30
  sp_rates[29] -= (fwd_rates[159] - rev_rates[145]);
  //sp 31
  sp_rates[30] += (fwd_rates[159] - rev_rates[145]);

  //rxn 160
  //sp 19
  sp_rates[18] -= (fwd_rates[160] - rev_rates[146]);
  //sp 30
  sp_rates[29] -= (fwd_rates[160] - rev_rates[146]);
  //sp 31
  sp_rates[30] += (fwd_rates[160] - rev_rates[146]);
  //sp 16
  sp_rates[15] += (fwd_rates[160] - rev_rates[146]);

  //rxn 161
  //sp 17
  sp_rates[16] -= (fwd_rates[161] - rev_rates[147]);
  //sp 30
  sp_rates[29] -= (fwd_rates[161] - rev_rates[147]);
  //sp 15
  sp_rates[14] += (fwd_rates[161] - rev_rates[147]);
  //sp 31
  sp_rates[30] += (fwd_rates[161] - rev_rates[147]);

  //rxn 162
  //sp 14
  sp_rates[13] += (fwd_rates[162] - rev_rates[148]);
  //sp 30
  sp_rates[29] -= (fwd_rates[162] - rev_rates[148]);
  //sp 15
  sp_rates[14] -= (fwd_rates[162] - rev_rates[148]);
  //sp 31
  sp_rates[30] += (fwd_rates[162] - rev_rates[148]);

  //rxn 163
  //sp 22
  sp_rates[21] += (fwd_rates[163] - rev_rates[149]);
  //sp 30
  sp_rates[29] -= (fwd_rates[163] - rev_rates[149]);
  //sp 23
  sp_rates[22] -= (fwd_rates[163] - rev_rates[149]);
  //sp 31
  sp_rates[30] += (fwd_rates[163] - rev_rates[149]);

  //rxn 164
  //sp 20
  sp_rates[19] += (fwd_rates[164] - rev_rates[150]);
  //sp 21
  sp_rates[20] -= (fwd_rates[164] - rev_rates[150]);
  //sp 30
  sp_rates[29] -= (fwd_rates[164] - rev_rates[150]);
  //sp 31
  sp_rates[30] += (fwd_rates[164] - rev_rates[150]);

  //rxn 165
  //sp 33
  sp_rates[32] = -(fwd_rates[165] - rev_rates[151]);
  //sp 30
  sp_rates[29] -= (fwd_rates[165] - rev_rates[151]);
  //sp 31
  sp_rates[30] += (fwd_rates[165] - rev_rates[151]);
  //sp 32
  sp_rates[31] = (fwd_rates[165] - rev_rates[151]);

  //rxn 166
  //sp 34
  sp_rates[33] = -(fwd_rates[166] - rev_rates[152]) * pres_mod[20];
  //sp 4
  sp_rates[3] -= (fwd_rates[166] - rev_rates[152]) * pres_mod[20];
  //sp 31
  sp_rates[30] += (fwd_rates[166] - rev_rates[152]) * pres_mod[20];

  //rxn 167
  //sp 34
  sp_rates[33] += (fwd_rates[167] - rev_rates[153]);
  //sp 3
  sp_rates[2] += (fwd_rates[167] - rev_rates[153]);
  //sp 4
  sp_rates[3] -= (fwd_rates[167] - rev_rates[153]);
  //sp 31
  sp_rates[30] -= (fwd_rates[167] - rev_rates[153]);

  //rxn 168
  //sp 34
  sp_rates[33] -= 2.0 * (fwd_rates[168] - rev_rates[154]);
  //sp 35
  sp_rates[34] = (fwd_rates[168] - rev_rates[154]);
  //sp 31
  sp_rates[30] += (fwd_rates[168] - rev_rates[154]);

  //rxn 169
  //sp 34
  sp_rates[33] += (fwd_rates[169] - rev_rates[155]);
  //sp 14
  sp_rates[13] += (fwd_rates[169] - rev_rates[155]);
  //sp 15
  sp_rates[14] -= (fwd_rates[169] - rev_rates[155]);
  //sp 31
  sp_rates[30] -= (fwd_rates[169] - rev_rates[155]);

  //rxn 170
  //sp 4
  sp_rates[3] += (fwd_rates[170] - rev_rates[156]);
  //sp 45
  sp_rates[44] = (fwd_rates[170] - rev_rates[156]);
  //sp 6
  sp_rates[5] -= (fwd_rates[170] - rev_rates[156]);
  //sp 31
  sp_rates[30] -= (fwd_rates[170] - rev_rates[156]);

  //rxn 171
  //sp 39
  sp_rates[38] = (fwd_rates[171] - rev_rates[157]);
  //sp 10
  sp_rates[9] -= (fwd_rates[171] - rev_rates[157]);
  //sp 31
  sp_rates[30] -= (fwd_rates[171] - rev_rates[157]);
  //sp 8
  sp_rates[7] += (fwd_rates[171] - rev_rates[157]);

  //rxn 172
  //sp 39
  sp_rates[38] += (fwd_rates[172] - rev_rates[158]);
  //sp 21
  sp_rates[20] -= (fwd_rates[172] - rev_rates[158]);
  //sp 31
  sp_rates[30] -= (fwd_rates[172] - rev_rates[158]);
  //sp 23
  sp_rates[22] += (fwd_rates[172] - rev_rates[158]);

  //rxn 173
  //sp 4
  sp_rates[3] += (fwd_rates[173] - rev_rates[159]);
  //sp 15
  sp_rates[14] -= 2.0 * (fwd_rates[173] - rev_rates[159]);
  //sp 31
  sp_rates[30] += (fwd_rates[173] - rev_rates[159]);

  //rxn 174
  //sp 33
  sp_rates[32] += (fwd_rates[174] - rev_rates[160]);
  //sp 5
  sp_rates[4] -= (fwd_rates[174] - rev_rates[160]);
  //sp 31
  sp_rates[30] -= (fwd_rates[174] - rev_rates[160]);

  //rxn 175
  //sp 42
  sp_rates[41] = (fwd_rates[175] - rev_rates[161]);
  //sp 5
  sp_rates[4] -= (fwd_rates[175] - rev_rates[161]);
  //sp 31
  sp_rates[30] -= (fwd_rates[175] - rev_rates[161]);

  //rxn 176
  //sp 10
  sp_rates[9] += (fwd_rates[176] - rev_rates[162]);
  //sp 5
  sp_rates[4] -= (fwd_rates[176] - rev_rates[162]);
  //sp 31
  sp_rates[30] -= (fwd_rates[176] - rev_rates[162]);
  //sp 34
  sp_rates[33] += (fwd_rates[176] - rev_rates[162]);

  //rxn 177
  //sp 43
  sp_rates[42] = (fwd_rates[177] - rev_rates[163]);
  //sp 5
  sp_rates[4] -= (fwd_rates[177] - rev_rates[163]);
  //sp 31
  sp_rates[30] -= (fwd_rates[177] - rev_rates[163]);
  //sp 8
  sp_rates[7] += (fwd_rates[177] - rev_rates[163]);

  //rxn 178
  //sp 5
  sp_rates[4] -= (fwd_rates[178] - rev_rates[164]);
  //sp 45
  sp_rates[44] += (fwd_rates[178] - rev_rates[164]);
  //sp 31
  sp_rates[30] -= (fwd_rates[178] - rev_rates[164]);
  //sp 8
  sp_rates[7] += (fwd_rates[178] - rev_rates[164]);

  //rxn 179
  //sp 33
  sp_rates[32] += (fwd_rates[179] - rev_rates[165]);
  //sp 42
  sp_rates[41] -= (fwd_rates[179] - rev_rates[165]);

  //rxn 180
  //sp 33
  sp_rates[32] -= (fwd_rates[180] - rev_rates[166]);
  //sp 45
  sp_rates[44] += (fwd_rates[180] - rev_rates[166]);
  //sp 8
  sp_rates[7] += (fwd_rates[180] - rev_rates[166]);

  //rxn 181
  //sp 33
  sp_rates[32] -= (fwd_rates[181] - rev_rates[167]);
  //sp 10
  sp_rates[9] += (fwd_rates[181] - rev_rates[167]);
  //sp 34
  sp_rates[33] += (fwd_rates[181] - rev_rates[167]);

  //rxn 182
  //sp 33
  sp_rates[32] -= (fwd_rates[182] - rev_rates[168]);
  //sp 43
  sp_rates[42] += (fwd_rates[182] - rev_rates[168]);
  //sp 8
  sp_rates[7] += (fwd_rates[182] - rev_rates[168]);

  //rxn 183
  //sp 42
  sp_rates[41] -= (fwd_rates[183] - rev_rates[169]);
  //sp 43
  sp_rates[42] += (fwd_rates[183] - rev_rates[169]);
  //sp 8
  sp_rates[7] += (fwd_rates[183] - rev_rates[169]);

  //rxn 184
  //sp 42
  sp_rates[41] -= (fwd_rates[184] - rev_rates[170]);
  //sp 10
  sp_rates[9] += (fwd_rates[184] - rev_rates[170]);
  //sp 34
  sp_rates[33] += (fwd_rates[184] - rev_rates[170]);

  //rxn 185
  //sp 42
  sp_rates[41] -= (fwd_rates[185] - rev_rates[171]);
  //sp 45
  sp_rates[44] += (fwd_rates[185] - rev_rates[171]);
  //sp 8
  sp_rates[7] += (fwd_rates[185] - rev_rates[171]);

  //rxn 186
  //sp 33
  sp_rates[32] -= (fwd_rates[186] - rev_rates[172]);
  //sp 3
  sp_rates[2] -= (fwd_rates[186] - rev_rates[172]);
  //sp 4
  sp_rates[3] += (fwd_rates[186] - rev_rates[172]);
  //sp 32
  sp_rates[31] += (fwd_rates[186] - rev_rates[172]);

  //rxn 187
  //sp 33
  sp_rates[32] -= (fwd_rates[187] - rev_rates[173]);
  //sp 10
  sp_rates[9] -= (fwd_rates[187] - rev_rates[173]);
  //sp 5
  sp_rates[4] += (fwd_rates[187] - rev_rates[173]);
  //sp 32
  sp_rates[31] += (fwd_rates[187] - rev_rates[173]);

  //rxn 188
  //sp 25
  sp_rates[24] -= (fwd_rates[188] - rev_rates[174]);
  //sp 26
  sp_rates[25] += (fwd_rates[188] - rev_rates[174]);
  //sp 33
  sp_rates[32] -= (fwd_rates[188] - rev_rates[174]);
  //sp 32
  sp_rates[31] += (fwd_rates[188] - rev_rates[174]);

  //rxn 189
  //sp 33
  sp_rates[32] -= (fwd_rates[189] - rev_rates[175]);
  //sp 14
  sp_rates[13] -= (fwd_rates[189] - rev_rates[175]);
  //sp 15
  sp_rates[14] += (fwd_rates[189] - rev_rates[175]);
  //sp 32
  sp_rates[31] += (fwd_rates[189] - rev_rates[175]);

  //rxn 190
  //sp 33
  sp_rates[32] -= (fwd_rates[190] - rev_rates[176]);
  //sp 32
  sp_rates[31] += (fwd_rates[190] - rev_rates[176]);
  //sp 22
  sp_rates[21] -= (fwd_rates[190] - rev_rates[176]);
  //sp 24
  sp_rates[23] += (fwd_rates[190] - rev_rates[176]);

  //rxn 191
  //sp 8
  sp_rates[7] += (fwd_rates[191] - rev_rates[177]);
  //sp 39
  sp_rates[38] += (fwd_rates[191] - rev_rates[177]);
  //sp 32
  sp_rates[31] -= (fwd_rates[191] - rev_rates[177]);

  //rxn 192
  //sp 26
  sp_rates[25] += (fwd_rates[192] - rev_rates[178]);
  //sp 43
  sp_rates[42] -= (fwd_rates[192] - rev_rates[178]);
  //sp 15
  sp_rates[14] += (fwd_rates[192] - rev_rates[178]);

  //rxn 193
  //sp 43
  sp_rates[42] -= (fwd_rates[193] - rev_rates[179]);
  //sp 45
  sp_rates[44] += (fwd_rates[193] - rev_rates[179]);

  //rxn 194
  //sp 3
  sp_rates[2] += (fwd_rates[194] - rev_rates[180]);
  //sp 44
  sp_rates[43] = (fwd_rates[194] - rev_rates[180]);
  //sp 43
  sp_rates[42] -= (fwd_rates[194] - rev_rates[180]);
  //sp 4
  sp_rates[3] -= (fwd_rates[194] - rev_rates[180]);

  //rxn 195
  //sp 43
  sp_rates[42] -= (fwd_rates[195] - rev_rates[181]);
  //sp 44
  sp_rates[43] += (fwd_rates[195] - rev_rates[181]);
  //sp 7
  sp_rates[6] += (fwd_rates[195] - rev_rates[181]);
  //sp 8
  sp_rates[7] -= (fwd_rates[195] - rev_rates[181]);

  //rxn 196
  //sp 9
  sp_rates[8] += (fwd_rates[196] - rev_rates[182]);
  //sp 10
  sp_rates[9] -= (fwd_rates[196] - rev_rates[182]);
  //sp 43
  sp_rates[42] -= (fwd_rates[196] - rev_rates[182]);
  //sp 44
  sp_rates[43] += (fwd_rates[196] - rev_rates[182]);

  //rxn 197
  //sp 43
  sp_rates[42] -= (fwd_rates[197] - rev_rates[183]);
  //sp 44
  sp_rates[43] += (fwd_rates[197] - rev_rates[183]);
  //sp 14
  sp_rates[13] += (fwd_rates[197] - rev_rates[183]);
  //sp 15
  sp_rates[14] -= (fwd_rates[197] - rev_rates[183]);

  //rxn 198
  //sp 43
  sp_rates[42] -= (fwd_rates[198] - rev_rates[184]);
  //sp 44
  sp_rates[43] += (fwd_rates[198] - rev_rates[184]);
  //sp 22
  sp_rates[21] += (fwd_rates[198] - rev_rates[184]);
  //sp 23
  sp_rates[22] -= (fwd_rates[198] - rev_rates[184]);

  //rxn 199
  //sp 44
  sp_rates[43] += (fwd_rates[199] - rev_rates[185]);
  //sp 43
  sp_rates[42] -= (fwd_rates[199] - rev_rates[185]);
  //sp 20
  sp_rates[19] += (fwd_rates[199] - rev_rates[185]);
  //sp 21
  sp_rates[20] -= (fwd_rates[199] - rev_rates[185]);

  //rxn 200
  //sp 33
  sp_rates[32] -= (fwd_rates[200] - rev_rates[186]);
  //sp 43
  sp_rates[42] -= (fwd_rates[200] - rev_rates[186]);
  //sp 44
  sp_rates[43] += (fwd_rates[200] - rev_rates[186]);
  //sp 32
  sp_rates[31] += (fwd_rates[200] - rev_rates[186]);

  //rxn 201
  //sp 44
  sp_rates[43] -= (fwd_rates[201] - rev_rates[187]);
  //sp 46
  sp_rates[45] = (fwd_rates[201] - rev_rates[187]);

  //rxn 202
  //sp 44
  sp_rates[43] -= (fwd_rates[202] - rev_rates[188]);
  //sp 47
  sp_rates[46] = (fwd_rates[202] - rev_rates[188]);

  //rxn 203
  //sp 34
  sp_rates[33] += (fwd_rates[203] - rev_rates[189]) * pres_mod[21];
  //sp 35
  sp_rates[34] -= (fwd_rates[203] - rev_rates[189]) * pres_mod[21];
  //sp 4
  sp_rates[3] -= (fwd_rates[203] - rev_rates[189]) * pres_mod[21];

  //rxn 204
  //sp 34
  sp_rates[33] -= (fwd_rates[204] - rev_rates[190]) * pres_mod[22];
  //sp 3
  sp_rates[2] += (fwd_rates[204] - rev_rates[190]) * pres_mod[22];
  //sp 36
  sp_rates[35] = (fwd_rates[204] - rev_rates[190]) * pres_mod[22];

  //rxn 205
  //sp 34
  sp_rates[33] -= (fwd_rates[205] - rev_rates[191]);
  //sp 35
  sp_rates[34] += (fwd_rates[205] - rev_rates[191]);
  //sp 5
  sp_rates[4] -= (fwd_rates[205] - rev_rates[191]);
  //sp 10
  sp_rates[9] += (fwd_rates[205] - rev_rates[191]);

  //rxn 206
  //sp 34
  sp_rates[33] -= (fwd_rates[206] - rev_rates[192]);
  //sp 3
  sp_rates[2] += (fwd_rates[206] - rev_rates[192]);
  //sp 4
  sp_rates[3] -= (fwd_rates[206] - rev_rates[192]);
  //sp 35
  sp_rates[34] += (fwd_rates[206] - rev_rates[192]);

  //rxn 207
  //sp 34
  sp_rates[33] -= (fwd_rates[207] - rev_rates[193]);
  //sp 35
  sp_rates[34] += (fwd_rates[207] - rev_rates[193]);
  //sp 7
  sp_rates[6] += (fwd_rates[207] - rev_rates[193]);
  //sp 8
  sp_rates[7] -= (fwd_rates[207] - rev_rates[193]);

  //rxn 208
  //sp 34
  sp_rates[33] -= (fwd_rates[208] - rev_rates[194]);
  //sp 35
  sp_rates[34] += (fwd_rates[208] - rev_rates[194]);
  //sp 22
  sp_rates[21] += (fwd_rates[208] - rev_rates[194]);
  //sp 23
  sp_rates[22] -= (fwd_rates[208] - rev_rates[194]);

  //rxn 209
  //sp 34
  sp_rates[33] -= (fwd_rates[209] - rev_rates[195]);
  //sp 35
  sp_rates[34] += (fwd_rates[209] - rev_rates[195]);
  //sp 20
  sp_rates[19] += (fwd_rates[209] - rev_rates[195]);
  //sp 21
  sp_rates[20] -= (fwd_rates[209] - rev_rates[195]);

  //rxn 210
  //sp 33
  sp_rates[32] -= (fwd_rates[210] - rev_rates[196]);
  //sp 34
  sp_rates[33] -= (fwd_rates[210] - rev_rates[196]);
  //sp 35
  sp_rates[34] += (fwd_rates[210] - rev_rates[196]);
  //sp 32
  sp_rates[31] += (fwd_rates[210] - rev_rates[196]);

  //rxn 211
  //sp 34
  sp_rates[33] -= (fwd_rates[211] - rev_rates[197]);
  //sp 35
  sp_rates[34] += (fwd_rates[211] - rev_rates[197]);
  //sp 14
  sp_rates[13] += (fwd_rates[211] - rev_rates[197]);
  //sp 15
  sp_rates[14] -= (fwd_rates[211] - rev_rates[197]);

  //rxn 212
  //sp 34
  sp_rates[33] -= (fwd_rates[212] - rev_rates[198]);
  //sp 35
  sp_rates[34] += (fwd_rates[212] - rev_rates[198]);
  //sp 14
  sp_rates[13] += (fwd_rates[212] - rev_rates[198]);
  //sp 15
  sp_rates[14] -= (fwd_rates[212] - rev_rates[198]);

  //rxn 213
  //sp 34
  sp_rates[33] -= (fwd_rates[213] - rev_rates[199]);
  //sp 6
  sp_rates[5] -= (fwd_rates[213] - rev_rates[199]);
  //sp 15
  sp_rates[14] += (fwd_rates[213] - rev_rates[199]);
  //sp 26
  sp_rates[25] += (fwd_rates[213] - rev_rates[199]);

  //rxn 214
  //sp 34
  sp_rates[33] -= (fwd_rates[214] - rev_rates[200]);
  //sp 4
  sp_rates[3] += (fwd_rates[214] - rev_rates[200]);
  //sp 6
  sp_rates[5] -= (fwd_rates[214] - rev_rates[200]);
  //sp 47
  sp_rates[46] += (fwd_rates[214] - rev_rates[200]);

  //rxn 215
  //sp 34
  sp_rates[33] += (fwd_rates[215] - rev_rates[201]);
  //sp 19
  sp_rates[18] -= (fwd_rates[215] - rev_rates[201]);
  //sp 4
  sp_rates[3] += (fwd_rates[215] - rev_rates[201]);
  //sp 14
  sp_rates[13] -= (fwd_rates[215] - rev_rates[201]);

  //rxn 216
  //sp 17
  sp_rates[16] -= (fwd_rates[216] - rev_rates[202]);
  //sp 34
  sp_rates[33] += (fwd_rates[216] - rev_rates[202]);
  //sp 4
  sp_rates[3] += (fwd_rates[216] - rev_rates[202]);
  //sp 15
  sp_rates[14] -= (fwd_rates[216] - rev_rates[202]);

  //rxn 217
  //sp 25
  sp_rates[24] += (fwd_rates[217] - rev_rates[203]);
  //sp 34
  sp_rates[33] -= (fwd_rates[217] - rev_rates[203]);
  //sp 15
  sp_rates[14] += (fwd_rates[217] - rev_rates[203]);
  //sp 8
  sp_rates[7] -= (fwd_rates[217] - rev_rates[203]);

  //rxn 218
  //sp 34
  sp_rates[33] -= (fwd_rates[218] - rev_rates[204]);
  //sp 4
  sp_rates[3] += (fwd_rates[218] - rev_rates[204]);
  //sp 45
  sp_rates[44] += (fwd_rates[218] - rev_rates[204]);
  //sp 8
  sp_rates[7] -= (fwd_rates[218] - rev_rates[204]);

  //rxn 219
  //sp 34
  sp_rates[33] -= (fwd_rates[219] - rev_rates[205]);
  //sp 40
  sp_rates[39] = (fwd_rates[219] - rev_rates[205]);
  //sp 8
  sp_rates[7] -= (fwd_rates[219] - rev_rates[205]);

  //rxn 220
  //sp 10
  sp_rates[9] -= (fwd_rates[220] - rev_rates[206]);
  //sp 43
  sp_rates[42] += (fwd_rates[220] - rev_rates[206]);
  //sp 8
  sp_rates[7] += (fwd_rates[220] - rev_rates[206]);
  //sp 34
  sp_rates[33] -= (fwd_rates[220] - rev_rates[206]);

  //rxn 221
  //sp 34
  sp_rates[33] -= (fwd_rates[221] - rev_rates[207]);
  //sp 43
  sp_rates[42] += (fwd_rates[221] - rev_rates[207]);
  //sp 21
  sp_rates[20] -= (fwd_rates[221] - rev_rates[207]);
  //sp 23
  sp_rates[22] += (fwd_rates[221] - rev_rates[207]);

  //rxn 222
  //sp 33
  sp_rates[32] -= (fwd_rates[222] - rev_rates[208]);
  //sp 34
  sp_rates[33] -= (fwd_rates[222] - rev_rates[208]);
  //sp 43
  sp_rates[42] += (fwd_rates[222] - rev_rates[208]);
  //sp 39
  sp_rates[38] += (fwd_rates[222] - rev_rates[208]);

  //rxn 223
  //sp 36
  sp_rates[35] -= (fwd_rates[223] - rev_rates[209]) * pres_mod[23];
  //sp 35
  sp_rates[34] += (fwd_rates[223] - rev_rates[209]) * pres_mod[23];
  //sp 4
  sp_rates[3] -= (fwd_rates[223] - rev_rates[209]) * pres_mod[23];

  //rxn 224
  //sp 35
  sp_rates[34] -= (fwd_rates[224] - rev_rates[210]);
  //sp 5
  sp_rates[4] -= (fwd_rates[224] - rev_rates[210]);
  //sp 6
  sp_rates[5] += (fwd_rates[224] - rev_rates[210]);
  //sp 47
  sp_rates[46] += (fwd_rates[224] - rev_rates[210]);

  //rxn 225
  //sp 35
  sp_rates[34] -= (fwd_rates[225] - rev_rates[211]);
  //sp 5
  sp_rates[4] -= (fwd_rates[225] - rev_rates[211]);
  //sp 6
  sp_rates[5] += (fwd_rates[225] - rev_rates[211]);
  //sp 47
  sp_rates[46] += (fwd_rates[225] - rev_rates[211]);

  //rxn 226
  //sp 10
  sp_rates[9] += (fwd_rates[226] - rev_rates[212]);
  //sp 35
  sp_rates[34] -= (fwd_rates[226] - rev_rates[212]);
  //sp 36
  sp_rates[35] += (fwd_rates[226] - rev_rates[212]);
  //sp 5
  sp_rates[4] -= (fwd_rates[226] - rev_rates[212]);

  //rxn 227
  //sp 10
  sp_rates[9] += (fwd_rates[227] - rev_rates[213]);
  //sp 35
  sp_rates[34] -= (fwd_rates[227] - rev_rates[213]);
  //sp 36
  sp_rates[35] += (fwd_rates[227] - rev_rates[213]);
  //sp 5
  sp_rates[4] -= (fwd_rates[227] - rev_rates[213]);

  //rxn 228
  //sp 35
  sp_rates[34] -= (fwd_rates[228] - rev_rates[214]);
  //sp 5
  sp_rates[4] -= (fwd_rates[228] - rev_rates[214]);
  //sp 48
  sp_rates[47] = (fwd_rates[228] - rev_rates[214]);
  //sp 8
  sp_rates[7] += (fwd_rates[228] - rev_rates[214]);

  //rxn 229
  //sp 35
  sp_rates[34] -= (fwd_rates[229] - rev_rates[215]);
  //sp 5
  sp_rates[4] -= (fwd_rates[229] - rev_rates[215]);
  //sp 48
  sp_rates[47] += (fwd_rates[229] - rev_rates[215]);
  //sp 8
  sp_rates[7] += (fwd_rates[229] - rev_rates[215]);

  //rxn 230
  //sp 25
  sp_rates[24] += (fwd_rates[230] - rev_rates[216]);
  //sp 26
  sp_rates[25] += (fwd_rates[230] - rev_rates[216]);
  //sp 35
  sp_rates[34] -= (fwd_rates[230] - rev_rates[216]);
  //sp 5
  sp_rates[4] -= (fwd_rates[230] - rev_rates[216]);

  //rxn 231
  //sp 25
  sp_rates[24] += (fwd_rates[231] - rev_rates[217]);
  //sp 26
  sp_rates[25] += (fwd_rates[231] - rev_rates[217]);
  //sp 35
  sp_rates[34] -= (fwd_rates[231] - rev_rates[217]);
  //sp 5
  sp_rates[4] -= (fwd_rates[231] - rev_rates[217]);

  //rxn 232
  //sp 35
  sp_rates[34] -= fwd_rates[232];
  //sp 4
  sp_rates[3] += fwd_rates[232];
  //sp 5
  sp_rates[4] -= fwd_rates[232];
  //sp 11
  sp_rates[10] += fwd_rates[232];
  //sp 25
  sp_rates[24] += fwd_rates[232];

  //rxn 233
  //sp 35
  sp_rates[34] -= fwd_rates[233];
  //sp 4
  sp_rates[3] += fwd_rates[233];
  //sp 5
  sp_rates[4] -= fwd_rates[233];
  //sp 11
  sp_rates[10] += fwd_rates[233];
  //sp 25
  sp_rates[24] += fwd_rates[233];

  //rxn 234
  //sp 35
  sp_rates[34] -= (fwd_rates[234] - rev_rates[218]);
  //sp 11
  sp_rates[10] += (fwd_rates[234] - rev_rates[218]);
  //sp 5
  sp_rates[4] -= (fwd_rates[234] - rev_rates[218]);
  //sp 23
  sp_rates[22] += (fwd_rates[234] - rev_rates[218]);

  //rxn 235
  //sp 35
  sp_rates[34] -= (fwd_rates[235] - rev_rates[219]);
  //sp 11
  sp_rates[10] += (fwd_rates[235] - rev_rates[219]);
  //sp 5
  sp_rates[4] -= (fwd_rates[235] - rev_rates[219]);
  //sp 23
  sp_rates[22] += (fwd_rates[235] - rev_rates[219]);

  //rxn 236
  //sp 35
  sp_rates[34] -= (fwd_rates[236] - rev_rates[220]);
  //sp 12
  sp_rates[11] += (fwd_rates[236] - rev_rates[220]);
  //sp 5
  sp_rates[4] -= (fwd_rates[236] - rev_rates[220]);
  //sp 15
  sp_rates[14] += (fwd_rates[236] - rev_rates[220]);

  //rxn 237
  //sp 35
  sp_rates[34] -= (fwd_rates[237] - rev_rates[221]);
  //sp 12
  sp_rates[11] += (fwd_rates[237] - rev_rates[221]);
  //sp 5
  sp_rates[4] -= (fwd_rates[237] - rev_rates[221]);
  //sp 15
  sp_rates[14] += (fwd_rates[237] - rev_rates[221]);

  //rxn 238
  //sp 3
  sp_rates[2] += (fwd_rates[238] - rev_rates[222]);
  //sp 36
  sp_rates[35] += (fwd_rates[238] - rev_rates[222]);
  //sp 35
  sp_rates[34] -= (fwd_rates[238] - rev_rates[222]);
  //sp 4
  sp_rates[3] -= (fwd_rates[238] - rev_rates[222]);

  //rxn 239
  //sp 35
  sp_rates[34] -= (fwd_rates[239] - rev_rates[223]);
  //sp 36
  sp_rates[35] += (fwd_rates[239] - rev_rates[223]);
  //sp 7
  sp_rates[6] += (fwd_rates[239] - rev_rates[223]);
  //sp 8
  sp_rates[7] -= (fwd_rates[239] - rev_rates[223]);

  //rxn 240
  //sp 35
  sp_rates[34] -= (fwd_rates[240] - rev_rates[224]);
  //sp 36
  sp_rates[35] += (fwd_rates[240] - rev_rates[224]);
  //sp 14
  sp_rates[13] += (fwd_rates[240] - rev_rates[224]);
  //sp 15
  sp_rates[14] -= (fwd_rates[240] - rev_rates[224]);

  //rxn 241
  //sp 34
  sp_rates[33] += (fwd_rates[241] - rev_rates[225]);
  //sp 35
  sp_rates[34] -= 2.0 * (fwd_rates[241] - rev_rates[225]);
  //sp 36
  sp_rates[35] += (fwd_rates[241] - rev_rates[225]);

  //rxn 242
  //sp 36
  sp_rates[35] += (fwd_rates[242] - rev_rates[226]) * pres_mod[24];
  //sp 4
  sp_rates[3] -= (fwd_rates[242] - rev_rates[226]) * pres_mod[24];
  //sp 37
  sp_rates[36] = -(fwd_rates[242] - rev_rates[226]) * pres_mod[24];

  //rxn 243
  //sp 11
  sp_rates[10] += (fwd_rates[243] - rev_rates[227]);
  //sp 36
  sp_rates[35] -= (fwd_rates[243] - rev_rates[227]);
  //sp 6
  sp_rates[5] -= (fwd_rates[243] - rev_rates[227]);
  //sp 16
  sp_rates[15] += (fwd_rates[243] - rev_rates[227]);

  //rxn 244
  //sp 4
  sp_rates[3] += (fwd_rates[244] - rev_rates[228]);
  //sp 49
  sp_rates[48] = (fwd_rates[244] - rev_rates[228]);
  //sp 36
  sp_rates[35] -= (fwd_rates[244] - rev_rates[228]);
  //sp 6
  sp_rates[5] -= (fwd_rates[244] - rev_rates[228]);

  //rxn 245
  //sp 10
  sp_rates[9] -= (fwd_rates[245] - rev_rates[229]);
  //sp 36
  sp_rates[35] -= (fwd_rates[245] - rev_rates[229]);
  //sp 48
  sp_rates[47] += (fwd_rates[245] - rev_rates[229]);
  //sp 8
  sp_rates[7] += (fwd_rates[245] - rev_rates[229]);

  //rxn 246
  //sp 26
  sp_rates[25] -= (fwd_rates[246] - rev_rates[230]);
  //sp 11
  sp_rates[10] += (fwd_rates[246] - rev_rates[230]);
  //sp 36
  sp_rates[35] -= (fwd_rates[246] - rev_rates[230]);
  //sp 35
  sp_rates[34] += (fwd_rates[246] - rev_rates[230]);

  //rxn 247
  //sp 4
  sp_rates[3] += (fwd_rates[247] - rev_rates[231]);
  //sp 36
  sp_rates[35] -= (fwd_rates[247] - rev_rates[231]);
  //sp 80
  sp_rates[79] = (fwd_rates[247] - rev_rates[231]);
  //sp 16
  sp_rates[15] -= (fwd_rates[247] - rev_rates[231]);

  //rxn 248
  //sp 17
  sp_rates[16] -= (fwd_rates[248] - rev_rates[232]);
  //sp 4
  sp_rates[3] += (fwd_rates[248] - rev_rates[232]);
  //sp 36
  sp_rates[35] -= (fwd_rates[248] - rev_rates[232]);
  //sp 80
  sp_rates[79] += (fwd_rates[248] - rev_rates[232]);

  //rxn 249
  //sp 49
  sp_rates[48] -= (fwd_rates[249] - rev_rates[233]);
  //sp 11
  sp_rates[10] += (fwd_rates[249] - rev_rates[233]);
  //sp 36
  sp_rates[35] -= (fwd_rates[249] - rev_rates[233]);
  //sp 80
  sp_rates[79] += (fwd_rates[249] - rev_rates[233]);

  //rxn 250
  //sp 36
  sp_rates[35] -= (fwd_rates[250] - rev_rates[234]);
  //sp 37
  sp_rates[36] += (fwd_rates[250] - rev_rates[234]);
  //sp 7
  sp_rates[6] += (fwd_rates[250] - rev_rates[234]);
  //sp 8
  sp_rates[7] -= (fwd_rates[250] - rev_rates[234]);

  //rxn 251
  //sp 4
  sp_rates[3] += (fwd_rates[251] - rev_rates[235]);
  //sp 36
  sp_rates[35] -= (fwd_rates[251] - rev_rates[235]);
  //sp 48
  sp_rates[47] += (fwd_rates[251] - rev_rates[235]);
  //sp 8
  sp_rates[7] -= (fwd_rates[251] - rev_rates[235]);

  //rxn 252
  //sp 11
  sp_rates[10] += (fwd_rates[252] - rev_rates[236]);
  //sp 36
  sp_rates[35] -= (fwd_rates[252] - rev_rates[236]);
  //sp 15
  sp_rates[14] += (fwd_rates[252] - rev_rates[236]);
  //sp 8
  sp_rates[7] -= (fwd_rates[252] - rev_rates[236]);

  //rxn 253
  //sp 26
  sp_rates[25] += (fwd_rates[253] - rev_rates[237]);
  //sp 11
  sp_rates[10] += (fwd_rates[253] - rev_rates[237]);
  //sp 5
  sp_rates[4] -= (fwd_rates[253] - rev_rates[237]);
  //sp 37
  sp_rates[36] -= (fwd_rates[253] - rev_rates[237]);

  //rxn 254
  //sp 11
  sp_rates[10] += (fwd_rates[254] - rev_rates[238]);
  //sp 19
  sp_rates[18] += (fwd_rates[254] - rev_rates[238]);
  //sp 37
  sp_rates[36] -= (fwd_rates[254] - rev_rates[238]);
  //sp 6
  sp_rates[5] -= (fwd_rates[254] - rev_rates[238]);

  //rxn 255
  //sp 36
  sp_rates[35] += (fwd_rates[255] - rev_rates[239]);
  //sp 3
  sp_rates[2] -= (fwd_rates[255] - rev_rates[239]);
  //sp 4
  sp_rates[3] += (fwd_rates[255] - rev_rates[239]);
  //sp 37
  sp_rates[36] -= (fwd_rates[255] - rev_rates[239]);

  //rxn 256
  //sp 49
  sp_rates[48] += (fwd_rates[256] - rev_rates[240]);
  //sp 4
  sp_rates[3] += (fwd_rates[256] - rev_rates[240]);
  //sp 37
  sp_rates[36] -= (fwd_rates[256] - rev_rates[240]);
  //sp 8
  sp_rates[7] -= (fwd_rates[256] - rev_rates[240]);

  //rxn 257
  //sp 26
  sp_rates[25] += (fwd_rates[257] - rev_rates[241]) * pres_mod[25];
  //sp 45
  sp_rates[44] -= (fwd_rates[257] - rev_rates[241]) * pres_mod[25];
  //sp 15
  sp_rates[14] += (fwd_rates[257] - rev_rates[241]) * pres_mod[25];

  //rxn 258
  //sp 11
  sp_rates[10] += (fwd_rates[258] - rev_rates[242]) * pres_mod[26];
  //sp 45
  sp_rates[44] -= (fwd_rates[258] - rev_rates[242]) * pres_mod[26];
  //sp 14
  sp_rates[13] += (fwd_rates[258] - rev_rates[242]) * pres_mod[26];

  //rxn 259
  //sp 10
  sp_rates[9] += (fwd_rates[259] - rev_rates[243]);
  //sp 5
  sp_rates[4] -= (fwd_rates[259] - rev_rates[243]);
  //sp 46
  sp_rates[45] += (fwd_rates[259] - rev_rates[243]);
  //sp 45
  sp_rates[44] -= (fwd_rates[259] - rev_rates[243]);

  //rxn 260
  //sp 46
  sp_rates[45] += (fwd_rates[260] - rev_rates[244]);
  //sp 45
  sp_rates[44] -= (fwd_rates[260] - rev_rates[244]);
  //sp 6
  sp_rates[5] -= (fwd_rates[260] - rev_rates[244]);
  //sp 8
  sp_rates[7] += (fwd_rates[260] - rev_rates[244]);

  //rxn 261
  //sp 45
  sp_rates[44] -= (fwd_rates[261] - rev_rates[245]);
  //sp 46
  sp_rates[45] += (fwd_rates[261] - rev_rates[245]);
  //sp 7
  sp_rates[6] += (fwd_rates[261] - rev_rates[245]);
  //sp 8
  sp_rates[7] -= (fwd_rates[261] - rev_rates[245]);

  //rxn 262
  //sp 9
  sp_rates[8] += (fwd_rates[262] - rev_rates[246]);
  //sp 10
  sp_rates[9] -= (fwd_rates[262] - rev_rates[246]);
  //sp 45
  sp_rates[44] -= (fwd_rates[262] - rev_rates[246]);
  //sp 46
  sp_rates[45] += (fwd_rates[262] - rev_rates[246]);

  //rxn 263
  //sp 20
  sp_rates[19] += (fwd_rates[263] - rev_rates[247]);
  //sp 21
  sp_rates[20] -= (fwd_rates[263] - rev_rates[247]);
  //sp 46
  sp_rates[45] += (fwd_rates[263] - rev_rates[247]);
  //sp 45
  sp_rates[44] -= (fwd_rates[263] - rev_rates[247]);

  //rxn 264
  //sp 3
  sp_rates[2] += (fwd_rates[264] - rev_rates[248]);
  //sp 4
  sp_rates[3] -= (fwd_rates[264] - rev_rates[248]);
  //sp 45
  sp_rates[44] -= (fwd_rates[264] - rev_rates[248]);
  //sp 47
  sp_rates[46] += (fwd_rates[264] - rev_rates[248]);

  //rxn 265
  //sp 47
  sp_rates[46] += (fwd_rates[265] - rev_rates[249]);
  //sp 45
  sp_rates[44] -= (fwd_rates[265] - rev_rates[249]);
  //sp 7
  sp_rates[6] += (fwd_rates[265] - rev_rates[249]);
  //sp 8
  sp_rates[7] -= (fwd_rates[265] - rev_rates[249]);

  //rxn 266
  //sp 28
  sp_rates[27] += (fwd_rates[266] - rev_rates[250]);
  //sp 45
  sp_rates[44] -= (fwd_rates[266] - rev_rates[250]);
  //sp 15
  sp_rates[14] += (fwd_rates[266] - rev_rates[250]);
  //sp 8
  sp_rates[7] -= (fwd_rates[266] - rev_rates[250]);

  //rxn 267
  //sp 11
  sp_rates[10] += (fwd_rates[267] - rev_rates[251]);
  //sp 46
  sp_rates[45] -= (fwd_rates[267] - rev_rates[251]);
  //sp 15
  sp_rates[14] += (fwd_rates[267] - rev_rates[251]);

  //rxn 268
  //sp 51
  sp_rates[50] = (fwd_rates[268] - rev_rates[252]);
  //sp 5
  sp_rates[4] -= (fwd_rates[268] - rev_rates[252]);
  //sp 46
  sp_rates[45] -= (fwd_rates[268] - rev_rates[252]);

  //rxn 269
  //sp 10
  sp_rates[9] += (fwd_rates[269] - rev_rates[253]);
  //sp 5
  sp_rates[4] -= (fwd_rates[269] - rev_rates[253]);
  //sp 46
  sp_rates[45] -= (fwd_rates[269] - rev_rates[253]);
  //sp 48
  sp_rates[47] += (fwd_rates[269] - rev_rates[253]);

  //rxn 270
  //sp 4
  sp_rates[3] += (fwd_rates[270] - rev_rates[254]) * pres_mod[27];
  //sp 46
  sp_rates[45] -= (fwd_rates[270] - rev_rates[254]) * pres_mod[27];
  //sp 48
  sp_rates[47] += (fwd_rates[270] - rev_rates[254]) * pres_mod[27];

  //rxn 271
  //sp 3
  sp_rates[2] += (fwd_rates[271] - rev_rates[255]);
  //sp 4
  sp_rates[3] -= (fwd_rates[271] - rev_rates[255]);
  //sp 46
  sp_rates[45] -= (fwd_rates[271] - rev_rates[255]);
  //sp 48
  sp_rates[47] += (fwd_rates[271] - rev_rates[255]);

  //rxn 272
  //sp 48
  sp_rates[47] += (fwd_rates[272] - rev_rates[256]);
  //sp 46
  sp_rates[45] -= (fwd_rates[272] - rev_rates[256]);
  //sp 6
  sp_rates[5] -= (fwd_rates[272] - rev_rates[256]);
  //sp 8
  sp_rates[7] += (fwd_rates[272] - rev_rates[256]);

  //rxn 273
  //sp 14
  sp_rates[13] += (fwd_rates[273] - rev_rates[257]);
  //sp 46
  sp_rates[45] -= (fwd_rates[273] - rev_rates[257]);
  //sp 15
  sp_rates[14] -= (fwd_rates[273] - rev_rates[257]);
  //sp 48
  sp_rates[47] += (fwd_rates[273] - rev_rates[257]);

  //rxn 274
  //sp 4
  sp_rates[3] += (fwd_rates[274] - rev_rates[258]) * pres_mod[28];
  //sp 47
  sp_rates[46] -= (fwd_rates[274] - rev_rates[258]) * pres_mod[28];
  //sp 48
  sp_rates[47] += (fwd_rates[274] - rev_rates[258]) * pres_mod[28];

  //rxn 275
  //sp 11
  sp_rates[10] += (fwd_rates[275] - rev_rates[259]) * pres_mod[29];
  //sp 47
  sp_rates[46] -= (fwd_rates[275] - rev_rates[259]) * pres_mod[29];
  //sp 15
  sp_rates[14] += (fwd_rates[275] - rev_rates[259]) * pres_mod[29];

  //rxn 276
  //sp 11
  sp_rates[10] += (fwd_rates[276] - rev_rates[260]) * pres_mod[30];
  //sp 16
  sp_rates[15] += (fwd_rates[276] - rev_rates[260]) * pres_mod[30];
  //sp 48
  sp_rates[47] -= (fwd_rates[276] - rev_rates[260]) * pres_mod[30];

  //rxn 277
  //sp 35
  sp_rates[34] -= (fwd_rates[277] - rev_rates[261]);
  //sp 4
  sp_rates[3] += (fwd_rates[277] - rev_rates[261]);
  //sp 6
  sp_rates[5] -= (fwd_rates[277] - rev_rates[261]);
  //sp 48
  sp_rates[47] += (fwd_rates[277] - rev_rates[261]);

  //rxn 278
  //sp 11
  sp_rates[10] += (fwd_rates[278] - rev_rates[262]);
  //sp 4
  sp_rates[3] -= (fwd_rates[278] - rev_rates[262]);
  //sp 15
  sp_rates[14] += (fwd_rates[278] - rev_rates[262]);
  //sp 48
  sp_rates[47] -= (fwd_rates[278] - rev_rates[262]);

  //rxn 279
  //sp 26
  sp_rates[25] -= (fwd_rates[279] - rev_rates[263]);
  //sp 3
  sp_rates[2] += (fwd_rates[279] - rev_rates[263]);
  //sp 15
  sp_rates[14] -= (fwd_rates[279] - rev_rates[263]);
  //sp 48
  sp_rates[47] += (fwd_rates[279] - rev_rates[263]);

  //rxn 280
  //sp 25
  sp_rates[24] -= (fwd_rates[280] - rev_rates[264]);
  //sp 19
  sp_rates[18] -= (fwd_rates[280] - rev_rates[264]);
  //sp 4
  sp_rates[3] += (fwd_rates[280] - rev_rates[264]);
  //sp 48
  sp_rates[47] += (fwd_rates[280] - rev_rates[264]);

  //rxn 281
  //sp 12
  sp_rates[11] += (fwd_rates[281] - rev_rates[265]);
  //sp 16
  sp_rates[15] += (fwd_rates[281] - rev_rates[265]);
  //sp 6
  sp_rates[5] -= (fwd_rates[281] - rev_rates[265]);
  //sp 48
  sp_rates[47] -= (fwd_rates[281] - rev_rates[265]);

  //rxn 282
  //sp 11
  sp_rates[10] += (fwd_rates[282] - rev_rates[266]);
  //sp 24
  sp_rates[23] += (fwd_rates[282] - rev_rates[266]);
  //sp 48
  sp_rates[47] -= (fwd_rates[282] - rev_rates[266]);
  //sp 8
  sp_rates[7] -= (fwd_rates[282] - rev_rates[266]);

  //rxn 283
  //sp 17
  sp_rates[16] -= (fwd_rates[283] - rev_rates[267]);
  //sp 34
  sp_rates[33] += (fwd_rates[283] - rev_rates[267]);
  //sp 11
  sp_rates[10] += (fwd_rates[283] - rev_rates[267]);
  //sp 48
  sp_rates[47] -= (fwd_rates[283] - rev_rates[267]);

  //rxn 284
  //sp 11
  sp_rates[10] += (fwd_rates[284] - rev_rates[268]);
  //sp 31
  sp_rates[30] += (fwd_rates[284] - rev_rates[268]);
  //sp 15
  sp_rates[14] -= (fwd_rates[284] - rev_rates[268]);
  //sp 48
  sp_rates[47] -= (fwd_rates[284] - rev_rates[268]);

  //rxn 285
  //sp 49
  sp_rates[48] -= fwd_rates[285];
  //sp 3
  sp_rates[2] += fwd_rates[285];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[285];
  //sp 8
  sp_rates[7] -= fwd_rates[285];

  //rxn 286
  //sp 49
  sp_rates[48] -= fwd_rates[286];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[286];
  //sp 4
  sp_rates[3] += fwd_rates[286];
  //sp 6
  sp_rates[5] -= fwd_rates[286];

  //rxn 287
  //sp 49
  sp_rates[48] += (fwd_rates[287] - rev_rates[269]) * pres_mod[31];
  //sp 11
  sp_rates[10] -= (fwd_rates[287] - rev_rates[269]) * pres_mod[31];
  //sp 19
  sp_rates[18] -= (fwd_rates[287] - rev_rates[269]) * pres_mod[31];

  //rxn 288
  //sp 49
  sp_rates[48] -= (fwd_rates[288] - rev_rates[270]);
  //sp 17
  sp_rates[16] += (fwd_rates[288] - rev_rates[270]);
  //sp 11
  sp_rates[10] += (fwd_rates[288] - rev_rates[270]);
  //sp 4
  sp_rates[3] -= (fwd_rates[288] - rev_rates[270]);

  //rxn 289
  //sp 49
  sp_rates[48] -= fwd_rates[289];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[289];
  //sp 5
  sp_rates[4] -= fwd_rates[289];
  //sp 8
  sp_rates[7] += fwd_rates[289];

  //rxn 290
  //sp 4
  sp_rates[3] += fwd_rates[290];
  //sp 5
  sp_rates[4] -= fwd_rates[290];
  //sp 11
  sp_rates[10] += fwd_rates[290];
  //sp 12
  sp_rates[11] += fwd_rates[290];
  //sp 49
  sp_rates[48] -= fwd_rates[290];

  //rxn 291
  //sp 49
  sp_rates[48] -= (fwd_rates[291] - rev_rates[271]);
  //sp 36
  sp_rates[35] += (fwd_rates[291] - rev_rates[271]);
  //sp 19
  sp_rates[18] -= (fwd_rates[291] - rev_rates[271]);
  //sp 11
  sp_rates[10] += (fwd_rates[291] - rev_rates[271]);

  //rxn 292
  //sp 34
  sp_rates[33] += (fwd_rates[292] - rev_rates[272]);
  //sp 38
  sp_rates[37] = -(fwd_rates[292] - rev_rates[272]);
  //sp 7
  sp_rates[6] += (fwd_rates[292] - rev_rates[272]);

  //rxn 293
  //sp 38
  sp_rates[37] -= (fwd_rates[293] - rev_rates[273]);
  //sp 15
  sp_rates[14] += (fwd_rates[293] - rev_rates[273]);
  //sp 24
  sp_rates[23] += (fwd_rates[293] - rev_rates[273]);

  //rxn 294
  //sp 38
  sp_rates[37] -= (fwd_rates[294] - rev_rates[274]);
  //sp 31
  sp_rates[30] += (fwd_rates[294] - rev_rates[274]);
  //sp 8
  sp_rates[7] += (fwd_rates[294] - rev_rates[274]);

  //rxn 295
  //sp 38
  sp_rates[37] -= (fwd_rates[295] - rev_rates[275]);
  //sp 6
  sp_rates[5] -= (fwd_rates[295] - rev_rates[275]);
  //sp 39
  sp_rates[38] += (fwd_rates[295] - rev_rates[275]);
  //sp 8
  sp_rates[7] += (fwd_rates[295] - rev_rates[275]);

  //rxn 296
  //sp 39
  sp_rates[38] += (fwd_rates[296] - rev_rates[276]);
  //sp 38
  sp_rates[37] -= (fwd_rates[296] - rev_rates[276]);
  //sp 7
  sp_rates[6] += (fwd_rates[296] - rev_rates[276]);
  //sp 8
  sp_rates[7] -= (fwd_rates[296] - rev_rates[276]);

  //rxn 297
  //sp 3
  sp_rates[2] += (fwd_rates[297] - rev_rates[277]);
  //sp 4
  sp_rates[3] -= (fwd_rates[297] - rev_rates[277]);
  //sp 38
  sp_rates[37] -= (fwd_rates[297] - rev_rates[277]);
  //sp 39
  sp_rates[38] += (fwd_rates[297] - rev_rates[277]);

  //rxn 298
  //sp 9
  sp_rates[8] += (fwd_rates[298] - rev_rates[278]);
  //sp 10
  sp_rates[9] -= (fwd_rates[298] - rev_rates[278]);
  //sp 38
  sp_rates[37] -= (fwd_rates[298] - rev_rates[278]);
  //sp 39
  sp_rates[38] += (fwd_rates[298] - rev_rates[278]);

  //rxn 299
  //sp 39
  sp_rates[38] += (fwd_rates[299] - rev_rates[279]);
  //sp 14
  sp_rates[13] += (fwd_rates[299] - rev_rates[279]);
  //sp 38
  sp_rates[37] -= (fwd_rates[299] - rev_rates[279]);
  //sp 15
  sp_rates[14] -= (fwd_rates[299] - rev_rates[279]);

  //rxn 300
  //sp 20
  sp_rates[19] += (fwd_rates[300] - rev_rates[280]);
  //sp 21
  sp_rates[20] -= (fwd_rates[300] - rev_rates[280]);
  //sp 38
  sp_rates[37] -= (fwd_rates[300] - rev_rates[280]);
  //sp 39
  sp_rates[38] += (fwd_rates[300] - rev_rates[280]);

  //rxn 301
  //sp 30
  sp_rates[29] += (fwd_rates[301] - rev_rates[281]);
  //sp 38
  sp_rates[37] -= (fwd_rates[301] - rev_rates[281]);
  //sp 31
  sp_rates[30] -= (fwd_rates[301] - rev_rates[281]);
  //sp 40
  sp_rates[39] += (fwd_rates[301] - rev_rates[281]);

  //rxn 302
  //sp 41
  sp_rates[40] = (fwd_rates[302] - rev_rates[282]);
  //sp 10
  sp_rates[9] += (fwd_rates[302] - rev_rates[282]);
  //sp 5
  sp_rates[4] -= (fwd_rates[302] - rev_rates[282]);
  //sp 38
  sp_rates[37] -= (fwd_rates[302] - rev_rates[282]);

  //rxn 303
  //sp 41
  sp_rates[40] += (fwd_rates[303] - rev_rates[283]);
  //sp 38
  sp_rates[37] -= (fwd_rates[303] - rev_rates[283]);
  //sp 6
  sp_rates[5] -= (fwd_rates[303] - rev_rates[283]);
  //sp 8
  sp_rates[7] += (fwd_rates[303] - rev_rates[283]);

  //rxn 304
  //sp 41
  sp_rates[40] += (fwd_rates[304] - rev_rates[284]);
  //sp 3
  sp_rates[2] += (fwd_rates[304] - rev_rates[284]);
  //sp 4
  sp_rates[3] -= (fwd_rates[304] - rev_rates[284]);
  //sp 38
  sp_rates[37] -= (fwd_rates[304] - rev_rates[284]);

  //rxn 305
  //sp 41
  sp_rates[40] += (fwd_rates[305] - rev_rates[285]);
  //sp 38
  sp_rates[37] -= (fwd_rates[305] - rev_rates[285]);
  //sp 7
  sp_rates[6] += (fwd_rates[305] - rev_rates[285]);
  //sp 8
  sp_rates[7] -= (fwd_rates[305] - rev_rates[285]);

  //rxn 306
  //sp 9
  sp_rates[8] += (fwd_rates[306] - rev_rates[286]);
  //sp 10
  sp_rates[9] -= (fwd_rates[306] - rev_rates[286]);
  //sp 41
  sp_rates[40] += (fwd_rates[306] - rev_rates[286]);
  //sp 38
  sp_rates[37] -= (fwd_rates[306] - rev_rates[286]);

  //rxn 307
  //sp 41
  sp_rates[40] += (fwd_rates[307] - rev_rates[287]);
  //sp 14
  sp_rates[13] += (fwd_rates[307] - rev_rates[287]);
  //sp 38
  sp_rates[37] -= (fwd_rates[307] - rev_rates[287]);
  //sp 15
  sp_rates[14] -= (fwd_rates[307] - rev_rates[287]);

  //rxn 308
  //sp 41
  sp_rates[40] += (fwd_rates[308] - rev_rates[288]);
  //sp 20
  sp_rates[19] += (fwd_rates[308] - rev_rates[288]);
  //sp 21
  sp_rates[20] -= (fwd_rates[308] - rev_rates[288]);
  //sp 38
  sp_rates[37] -= (fwd_rates[308] - rev_rates[288]);

  //rxn 309
  //sp 41
  sp_rates[40] += (fwd_rates[309] - rev_rates[289]);
  //sp 30
  sp_rates[29] += (fwd_rates[309] - rev_rates[289]);
  //sp 38
  sp_rates[37] -= (fwd_rates[309] - rev_rates[289]);
  //sp 31
  sp_rates[30] -= (fwd_rates[309] - rev_rates[289]);

  //rxn 310
  //sp 10
  sp_rates[9] += (fwd_rates[310] - rev_rates[290]);
  //sp 5
  sp_rates[4] -= (fwd_rates[310] - rev_rates[290]);
  //sp 38
  sp_rates[37] -= (fwd_rates[310] - rev_rates[290]);
  //sp 40
  sp_rates[39] += (fwd_rates[310] - rev_rates[290]);

  //rxn 311
  //sp 40
  sp_rates[39] += (fwd_rates[311] - rev_rates[291]);
  //sp 38
  sp_rates[37] -= (fwd_rates[311] - rev_rates[291]);
  //sp 6
  sp_rates[5] -= (fwd_rates[311] - rev_rates[291]);
  //sp 8
  sp_rates[7] += (fwd_rates[311] - rev_rates[291]);

  //rxn 312
  //sp 3
  sp_rates[2] += (fwd_rates[312] - rev_rates[292]);
  //sp 4
  sp_rates[3] -= (fwd_rates[312] - rev_rates[292]);
  //sp 38
  sp_rates[37] -= (fwd_rates[312] - rev_rates[292]);
  //sp 40
  sp_rates[39] += (fwd_rates[312] - rev_rates[292]);

  //rxn 313
  //sp 40
  sp_rates[39] += (fwd_rates[313] - rev_rates[293]);
  //sp 38
  sp_rates[37] -= (fwd_rates[313] - rev_rates[293]);
  //sp 7
  sp_rates[6] += (fwd_rates[313] - rev_rates[293]);
  //sp 8
  sp_rates[7] -= (fwd_rates[313] - rev_rates[293]);

  //rxn 314
  //sp 9
  sp_rates[8] += (fwd_rates[314] - rev_rates[294]);
  //sp 10
  sp_rates[9] -= (fwd_rates[314] - rev_rates[294]);
  //sp 38
  sp_rates[37] -= (fwd_rates[314] - rev_rates[294]);
  //sp 40
  sp_rates[39] += (fwd_rates[314] - rev_rates[294]);

  //rxn 315
  //sp 14
  sp_rates[13] += (fwd_rates[315] - rev_rates[295]);
  //sp 38
  sp_rates[37] -= (fwd_rates[315] - rev_rates[295]);
  //sp 15
  sp_rates[14] -= (fwd_rates[315] - rev_rates[295]);
  //sp 40
  sp_rates[39] += (fwd_rates[315] - rev_rates[295]);

  //rxn 316
  //sp 20
  sp_rates[19] += (fwd_rates[316] - rev_rates[296]);
  //sp 21
  sp_rates[20] -= (fwd_rates[316] - rev_rates[296]);
  //sp 38
  sp_rates[37] -= (fwd_rates[316] - rev_rates[296]);
  //sp 40
  sp_rates[39] += (fwd_rates[316] - rev_rates[296]);

  //rxn 317
  //sp 41
  sp_rates[40] -= (fwd_rates[317] - rev_rates[297]);
  //sp 4
  sp_rates[3] += (fwd_rates[317] - rev_rates[297]);
  //sp 45
  sp_rates[44] += (fwd_rates[317] - rev_rates[297]);

  //rxn 318
  //sp 41
  sp_rates[40] -= (fwd_rates[318] - rev_rates[298]);
  //sp 39
  sp_rates[38] += (fwd_rates[318] - rev_rates[298]);

  //rxn 319
  //sp 41
  sp_rates[40] -= (fwd_rates[319] - rev_rates[299]);
  //sp 40
  sp_rates[39] += (fwd_rates[319] - rev_rates[299]);

  //rxn 320
  //sp 25
  sp_rates[24] -= (fwd_rates[320] - rev_rates[300]);
  //sp 39
  sp_rates[38] += (fwd_rates[320] - rev_rates[300]);
  //sp 15
  sp_rates[14] -= (fwd_rates[320] - rev_rates[300]);

  //rxn 321
  //sp 4
  sp_rates[3] -= (fwd_rates[321] - rev_rates[301]);
  //sp 45
  sp_rates[44] -= (fwd_rates[321] - rev_rates[301]);
  //sp 39
  sp_rates[38] += (fwd_rates[321] - rev_rates[301]);

  //rxn 322
  //sp 10
  sp_rates[9] += (fwd_rates[322] - rev_rates[302]);
  //sp 5
  sp_rates[4] -= (fwd_rates[322] - rev_rates[302]);
  //sp 39
  sp_rates[38] -= (fwd_rates[322] - rev_rates[302]);
  //sp 45
  sp_rates[44] += (fwd_rates[322] - rev_rates[302]);

  //rxn 323
  //sp 41
  sp_rates[40] -= (fwd_rates[323] - rev_rates[303]);
  //sp 10
  sp_rates[9] += (fwd_rates[323] - rev_rates[303]);
  //sp 5
  sp_rates[4] -= (fwd_rates[323] - rev_rates[303]);
  //sp 45
  sp_rates[44] += (fwd_rates[323] - rev_rates[303]);

  //rxn 324
  //sp 59
  sp_rates[58] = -(fwd_rates[324] - rev_rates[304]) * pres_mod[32];
  //sp 15
  sp_rates[14] += (fwd_rates[324] - rev_rates[304]) * pres_mod[32];
  //sp 31
  sp_rates[30] += (fwd_rates[324] - rev_rates[304]) * pres_mod[32];

  //rxn 325
  //sp 59
  sp_rates[58] += (fwd_rates[325] - rev_rates[305]);
  //sp 4
  sp_rates[3] -= (fwd_rates[325] - rev_rates[305]);
  //sp 61
  sp_rates[60] = -(fwd_rates[325] - rev_rates[305]);

  //rxn 326
  //sp 60
  sp_rates[59] = -(fwd_rates[326] - rev_rates[306]);
  //sp 59
  sp_rates[58] += (fwd_rates[326] - rev_rates[306]);
  //sp 4
  sp_rates[3] -= (fwd_rates[326] - rev_rates[306]);

  //rxn 327
  //sp 60
  sp_rates[59] -= (fwd_rates[327] - rev_rates[307]);
  //sp 61
  sp_rates[60] += (fwd_rates[327] - rev_rates[307]);

  //rxn 328
  //sp 10
  sp_rates[9] += (fwd_rates[328] - rev_rates[308]);
  //sp 59
  sp_rates[58] -= (fwd_rates[328] - rev_rates[308]);
  //sp 60
  sp_rates[59] += (fwd_rates[328] - rev_rates[308]);
  //sp 5
  sp_rates[4] -= (fwd_rates[328] - rev_rates[308]);

  //rxn 329
  //sp 3
  sp_rates[2] += (fwd_rates[329] - rev_rates[309]);
  //sp 60
  sp_rates[59] += (fwd_rates[329] - rev_rates[309]);
  //sp 59
  sp_rates[58] -= (fwd_rates[329] - rev_rates[309]);
  //sp 4
  sp_rates[3] -= (fwd_rates[329] - rev_rates[309]);

  //rxn 330
  //sp 3
  sp_rates[2] += (fwd_rates[330] - rev_rates[310]);
  //sp 59
  sp_rates[58] -= (fwd_rates[330] - rev_rates[310]);
  //sp 4
  sp_rates[3] -= (fwd_rates[330] - rev_rates[310]);
  //sp 61
  sp_rates[60] += (fwd_rates[330] - rev_rates[310]);

  //rxn 331
  //sp 59
  sp_rates[58] -= (fwd_rates[331] - rev_rates[311]);
  //sp 60
  sp_rates[59] += (fwd_rates[331] - rev_rates[311]);
  //sp 6
  sp_rates[5] -= (fwd_rates[331] - rev_rates[311]);
  //sp 8
  sp_rates[7] += (fwd_rates[331] - rev_rates[311]);

  //rxn 332
  //sp 59
  sp_rates[58] -= (fwd_rates[332] - rev_rates[312]);
  //sp 60
  sp_rates[59] += (fwd_rates[332] - rev_rates[312]);
  //sp 7
  sp_rates[6] += (fwd_rates[332] - rev_rates[312]);
  //sp 8
  sp_rates[7] -= (fwd_rates[332] - rev_rates[312]);

  //rxn 333
  //sp 9
  sp_rates[8] += (fwd_rates[333] - rev_rates[313]);
  //sp 10
  sp_rates[9] -= (fwd_rates[333] - rev_rates[313]);
  //sp 59
  sp_rates[58] -= (fwd_rates[333] - rev_rates[313]);
  //sp 60
  sp_rates[59] += (fwd_rates[333] - rev_rates[313]);

  //rxn 334
  //sp 59
  sp_rates[58] -= (fwd_rates[334] - rev_rates[314]);
  //sp 60
  sp_rates[59] += (fwd_rates[334] - rev_rates[314]);
  //sp 14
  sp_rates[13] += (fwd_rates[334] - rev_rates[314]);
  //sp 15
  sp_rates[14] -= (fwd_rates[334] - rev_rates[314]);

  //rxn 335
  //sp 34
  sp_rates[33] += (fwd_rates[335] - rev_rates[315]);
  //sp 35
  sp_rates[34] -= (fwd_rates[335] - rev_rates[315]);
  //sp 60
  sp_rates[59] += (fwd_rates[335] - rev_rates[315]);
  //sp 59
  sp_rates[58] -= (fwd_rates[335] - rev_rates[315]);

  //rxn 336
  //sp 59
  sp_rates[58] -= (fwd_rates[336] - rev_rates[316]);
  //sp 60
  sp_rates[59] += (fwd_rates[336] - rev_rates[316]);
  //sp 30
  sp_rates[29] += (fwd_rates[336] - rev_rates[316]);
  //sp 31
  sp_rates[30] -= (fwd_rates[336] - rev_rates[316]);

  //rxn 337
  //sp 59
  sp_rates[58] -= (fwd_rates[337] - rev_rates[317]);
  //sp 60
  sp_rates[59] += (fwd_rates[337] - rev_rates[317]);
  //sp 62
  sp_rates[61] = (fwd_rates[337] - rev_rates[317]);
  //sp 63
  sp_rates[62] = -(fwd_rates[337] - rev_rates[317]);

  //rxn 338
  //sp 59
  sp_rates[58] -= (fwd_rates[338] - rev_rates[318]);
  //sp 60
  sp_rates[59] += (fwd_rates[338] - rev_rates[318]);
  //sp 22
  sp_rates[21] += (fwd_rates[338] - rev_rates[318]);
  //sp 23
  sp_rates[22] -= (fwd_rates[338] - rev_rates[318]);

  //rxn 339
  //sp 60
  sp_rates[59] += (fwd_rates[339] - rev_rates[319]);
  //sp 59
  sp_rates[58] -= (fwd_rates[339] - rev_rates[319]);
  //sp 20
  sp_rates[19] += (fwd_rates[339] - rev_rates[319]);
  //sp 21
  sp_rates[20] -= (fwd_rates[339] - rev_rates[319]);

  //rxn 340
  //sp 33
  sp_rates[32] -= (fwd_rates[340] - rev_rates[320]);
  //sp 59
  sp_rates[58] -= (fwd_rates[340] - rev_rates[320]);
  //sp 60
  sp_rates[59] += (fwd_rates[340] - rev_rates[320]);
  //sp 32
  sp_rates[31] += (fwd_rates[340] - rev_rates[320]);

  //rxn 341
  //sp 10
  sp_rates[9] += (fwd_rates[341] - rev_rates[321]);
  //sp 59
  sp_rates[58] -= (fwd_rates[341] - rev_rates[321]);
  //sp 5
  sp_rates[4] -= (fwd_rates[341] - rev_rates[321]);
  //sp 61
  sp_rates[60] += (fwd_rates[341] - rev_rates[321]);

  //rxn 342
  //sp 59
  sp_rates[58] -= (fwd_rates[342] - rev_rates[322]);
  //sp 61
  sp_rates[60] += (fwd_rates[342] - rev_rates[322]);
  //sp 6
  sp_rates[5] -= (fwd_rates[342] - rev_rates[322]);
  //sp 8
  sp_rates[7] += (fwd_rates[342] - rev_rates[322]);

  //rxn 343
  //sp 59
  sp_rates[58] -= (fwd_rates[343] - rev_rates[323]);
  //sp 61
  sp_rates[60] += (fwd_rates[343] - rev_rates[323]);
  //sp 7
  sp_rates[6] += (fwd_rates[343] - rev_rates[323]);
  //sp 8
  sp_rates[7] -= (fwd_rates[343] - rev_rates[323]);

  //rxn 344
  //sp 9
  sp_rates[8] += (fwd_rates[344] - rev_rates[324]);
  //sp 10
  sp_rates[9] -= (fwd_rates[344] - rev_rates[324]);
  //sp 59
  sp_rates[58] -= (fwd_rates[344] - rev_rates[324]);
  //sp 61
  sp_rates[60] += (fwd_rates[344] - rev_rates[324]);

  //rxn 345
  //sp 59
  sp_rates[58] -= (fwd_rates[345] - rev_rates[325]);
  //sp 61
  sp_rates[60] += (fwd_rates[345] - rev_rates[325]);
  //sp 14
  sp_rates[13] += (fwd_rates[345] - rev_rates[325]);
  //sp 15
  sp_rates[14] -= (fwd_rates[345] - rev_rates[325]);

  //rxn 346
  //sp 34
  sp_rates[33] += (fwd_rates[346] - rev_rates[326]);
  //sp 35
  sp_rates[34] -= (fwd_rates[346] - rev_rates[326]);
  //sp 59
  sp_rates[58] -= (fwd_rates[346] - rev_rates[326]);
  //sp 61
  sp_rates[60] += (fwd_rates[346] - rev_rates[326]);

  //rxn 347
  //sp 59
  sp_rates[58] -= (fwd_rates[347] - rev_rates[327]);
  //sp 61
  sp_rates[60] += (fwd_rates[347] - rev_rates[327]);
  //sp 30
  sp_rates[29] += (fwd_rates[347] - rev_rates[327]);
  //sp 31
  sp_rates[30] -= (fwd_rates[347] - rev_rates[327]);

  //rxn 348
  //sp 59
  sp_rates[58] -= (fwd_rates[348] - rev_rates[328]);
  //sp 61
  sp_rates[60] += (fwd_rates[348] - rev_rates[328]);
  //sp 62
  sp_rates[61] += (fwd_rates[348] - rev_rates[328]);
  //sp 63
  sp_rates[62] -= (fwd_rates[348] - rev_rates[328]);

  //rxn 349
  //sp 59
  sp_rates[58] -= (fwd_rates[349] - rev_rates[329]);
  //sp 61
  sp_rates[60] += (fwd_rates[349] - rev_rates[329]);
  //sp 22
  sp_rates[21] += (fwd_rates[349] - rev_rates[329]);
  //sp 23
  sp_rates[22] -= (fwd_rates[349] - rev_rates[329]);

  //rxn 350
  //sp 59
  sp_rates[58] -= (fwd_rates[350] - rev_rates[330]);
  //sp 20
  sp_rates[19] += (fwd_rates[350] - rev_rates[330]);
  //sp 21
  sp_rates[20] -= (fwd_rates[350] - rev_rates[330]);
  //sp 61
  sp_rates[60] += (fwd_rates[350] - rev_rates[330]);

  //rxn 351
  //sp 33
  sp_rates[32] -= (fwd_rates[351] - rev_rates[331]);
  //sp 59
  sp_rates[58] -= (fwd_rates[351] - rev_rates[331]);
  //sp 61
  sp_rates[60] += (fwd_rates[351] - rev_rates[331]);
  //sp 32
  sp_rates[31] += (fwd_rates[351] - rev_rates[331]);

  //rxn 352
  //sp 60
  sp_rates[59] -= (fwd_rates[352] - rev_rates[332]);
  //sp 4
  sp_rates[3] -= (fwd_rates[352] - rev_rates[332]);
  //sp 15
  sp_rates[14] += (fwd_rates[352] - rev_rates[332]);
  //sp 31
  sp_rates[30] += (fwd_rates[352] - rev_rates[332]);

  //rxn 353
  //sp 60
  sp_rates[59] -= (fwd_rates[353] - rev_rates[333]);
  //sp 62
  sp_rates[61] += (fwd_rates[353] - rev_rates[333]);
  //sp 7
  sp_rates[6] += (fwd_rates[353] - rev_rates[333]);
  //sp 8
  sp_rates[7] -= (fwd_rates[353] - rev_rates[333]);

  //rxn 354
  //sp 4
  sp_rates[3] += (fwd_rates[354] - rev_rates[334]);
  //sp 60
  sp_rates[59] -= (fwd_rates[354] - rev_rates[334]);
  //sp 86
  sp_rates[85] = (fwd_rates[354] - rev_rates[334]);
  //sp 6
  sp_rates[5] -= (fwd_rates[354] - rev_rates[334]);

  //rxn 355
  //sp 60
  sp_rates[59] -= (fwd_rates[355] - rev_rates[335]);
  //sp 45
  sp_rates[44] += (fwd_rates[355] - rev_rates[335]);
  //sp 6
  sp_rates[5] -= (fwd_rates[355] - rev_rates[335]);
  //sp 15
  sp_rates[14] += (fwd_rates[355] - rev_rates[335]);

  //rxn 356
  //sp 8
  sp_rates[7] += fwd_rates[356];
  //sp 10
  sp_rates[9] -= fwd_rates[356];
  //sp 25
  sp_rates[24] += fwd_rates[356];
  //sp 61
  sp_rates[60] -= fwd_rates[356];
  //sp 31
  sp_rates[30] += fwd_rates[356];

  //rxn 357
  //sp 21
  sp_rates[20] -= fwd_rates[357];
  //sp 23
  sp_rates[22] += fwd_rates[357];
  //sp 25
  sp_rates[24] += fwd_rates[357];
  //sp 61
  sp_rates[60] -= fwd_rates[357];
  //sp 31
  sp_rates[30] += fwd_rates[357];

  //rxn 358
  //sp 8
  sp_rates[7] += fwd_rates[358];
  //sp 10
  sp_rates[9] -= fwd_rates[358];
  //sp 45
  sp_rates[44] += fwd_rates[358];
  //sp 15
  sp_rates[14] += fwd_rates[358];
  //sp 60
  sp_rates[59] -= fwd_rates[358];

  //rxn 359
  //sp 45
  sp_rates[44] += fwd_rates[359];
  //sp 15
  sp_rates[14] += fwd_rates[359];
  //sp 21
  sp_rates[20] -= fwd_rates[359];
  //sp 23
  sp_rates[22] += fwd_rates[359];
  //sp 60
  sp_rates[59] -= fwd_rates[359];

  //rxn 360
  //sp 35
  sp_rates[34] -= (fwd_rates[360] - rev_rates[336]) * pres_mod[33];
  //sp 62
  sp_rates[61] += (fwd_rates[360] - rev_rates[336]) * pres_mod[33];
  //sp 15
  sp_rates[14] -= (fwd_rates[360] - rev_rates[336]) * pres_mod[33];

  //rxn 361
  //sp 17
  sp_rates[16] -= (fwd_rates[361] - rev_rates[337]);
  //sp 34
  sp_rates[33] -= (fwd_rates[361] - rev_rates[337]);
  //sp 62
  sp_rates[61] += (fwd_rates[361] - rev_rates[337]);

  //rxn 362
  //sp 17
  sp_rates[16] -= (fwd_rates[362] - rev_rates[338]);
  //sp 34
  sp_rates[33] -= (fwd_rates[362] - rev_rates[338]);
  //sp 62
  sp_rates[61] += (fwd_rates[362] - rev_rates[338]);

  //rxn 363
  //sp 17
  sp_rates[16] -= (fwd_rates[363] - rev_rates[339]);
  //sp 34
  sp_rates[33] -= (fwd_rates[363] - rev_rates[339]);
  //sp 4
  sp_rates[3] += (fwd_rates[363] - rev_rates[339]);
  //sp 63
  sp_rates[62] += (fwd_rates[363] - rev_rates[339]);

  //rxn 364
  //sp 17
  sp_rates[16] -= (fwd_rates[364] - rev_rates[340]);
  //sp 34
  sp_rates[33] -= (fwd_rates[364] - rev_rates[340]);
  //sp 4
  sp_rates[3] += (fwd_rates[364] - rev_rates[340]);
  //sp 63
  sp_rates[62] += (fwd_rates[364] - rev_rates[340]);

  //rxn 365
  //sp 17
  sp_rates[16] -= (fwd_rates[365] - rev_rates[341]);
  //sp 34
  sp_rates[33] -= (fwd_rates[365] - rev_rates[341]);
  //sp 35
  sp_rates[34] += (fwd_rates[365] - rev_rates[341]);
  //sp 15
  sp_rates[14] += (fwd_rates[365] - rev_rates[341]);

  //rxn 366
  //sp 17
  sp_rates[16] -= (fwd_rates[366] - rev_rates[342]);
  //sp 34
  sp_rates[33] -= (fwd_rates[366] - rev_rates[342]);
  //sp 35
  sp_rates[34] += (fwd_rates[366] - rev_rates[342]);
  //sp 15
  sp_rates[14] += (fwd_rates[366] - rev_rates[342]);

  //rxn 367
  //sp 63
  sp_rates[62] += (fwd_rates[367] - rev_rates[343]);
  //sp 35
  sp_rates[34] -= (fwd_rates[367] - rev_rates[343]);
  //sp 4
  sp_rates[3] += (fwd_rates[367] - rev_rates[343]);
  //sp 15
  sp_rates[14] -= (fwd_rates[367] - rev_rates[343]);

  //rxn 368
  //sp 63
  sp_rates[62] += (fwd_rates[368] - rev_rates[344]);
  //sp 35
  sp_rates[34] -= (fwd_rates[368] - rev_rates[344]);
  //sp 4
  sp_rates[3] += (fwd_rates[368] - rev_rates[344]);
  //sp 15
  sp_rates[14] -= (fwd_rates[368] - rev_rates[344]);

  //rxn 369
  //sp 35
  sp_rates[34] += (fwd_rates[369] - rev_rates[345]);
  //sp 62
  sp_rates[61] -= (fwd_rates[369] - rev_rates[345]);
  //sp 15
  sp_rates[14] += (fwd_rates[369] - rev_rates[345]);

  //rxn 370
  //sp 35
  sp_rates[34] += (fwd_rates[370] - rev_rates[346]);
  //sp 62
  sp_rates[61] -= (fwd_rates[370] - rev_rates[346]);
  //sp 15
  sp_rates[14] += (fwd_rates[370] - rev_rates[346]);

  //rxn 371
  //sp 4
  sp_rates[3] += (fwd_rates[371] - rev_rates[347]);
  //sp 62
  sp_rates[61] -= (fwd_rates[371] - rev_rates[347]);
  //sp 63
  sp_rates[62] += (fwd_rates[371] - rev_rates[347]);

  //rxn 372
  //sp 4
  sp_rates[3] += (fwd_rates[372] - rev_rates[348]);
  //sp 62
  sp_rates[61] -= (fwd_rates[372] - rev_rates[348]);
  //sp 63
  sp_rates[62] += (fwd_rates[372] - rev_rates[348]);

  //rxn 373
  //sp 65
  sp_rates[64] = -(fwd_rates[373] - rev_rates[349]);
  //sp 4
  sp_rates[3] -= (fwd_rates[373] - rev_rates[349]);
  //sp 62
  sp_rates[61] += (fwd_rates[373] - rev_rates[349]);

  //rxn 374
  //sp 65
  sp_rates[64] -= (fwd_rates[374] - rev_rates[350]);
  //sp 4
  sp_rates[3] -= (fwd_rates[374] - rev_rates[350]);
  //sp 62
  sp_rates[61] += (fwd_rates[374] - rev_rates[350]);

  //rxn 375
  //sp 65
  sp_rates[64] -= (fwd_rates[375] - rev_rates[351]);
  //sp 63
  sp_rates[62] += (fwd_rates[375] - rev_rates[351]);

  //rxn 376
  //sp 65
  sp_rates[64] -= (fwd_rates[376] - rev_rates[352]);
  //sp 63
  sp_rates[62] += (fwd_rates[376] - rev_rates[352]);

  //rxn 377
  //sp 65
  sp_rates[64] -= (fwd_rates[377] - rev_rates[353]);
  //sp 35
  sp_rates[34] += (fwd_rates[377] - rev_rates[353]);
  //sp 4
  sp_rates[3] -= (fwd_rates[377] - rev_rates[353]);
  //sp 15
  sp_rates[14] += (fwd_rates[377] - rev_rates[353]);

  //rxn 378
  //sp 65
  sp_rates[64] -= (fwd_rates[378] - rev_rates[354]);
  //sp 35
  sp_rates[34] += (fwd_rates[378] - rev_rates[354]);
  //sp 4
  sp_rates[3] -= (fwd_rates[378] - rev_rates[354]);
  //sp 15
  sp_rates[14] += (fwd_rates[378] - rev_rates[354]);

  //rxn 379
  //sp 4
  sp_rates[3] += (fwd_rates[379] - rev_rates[355]);
  //sp 62
  sp_rates[61] -= (fwd_rates[379] - rev_rates[355]);
  //sp 64
  sp_rates[63] = (fwd_rates[379] - rev_rates[355]);

  //rxn 380
  //sp 3
  sp_rates[2] += (fwd_rates[380] - rev_rates[356]);
  //sp 4
  sp_rates[3] -= (fwd_rates[380] - rev_rates[356]);
  //sp 62
  sp_rates[61] -= (fwd_rates[380] - rev_rates[356]);
  //sp 63
  sp_rates[62] += (fwd_rates[380] - rev_rates[356]);

  //rxn 381
  //sp 10
  sp_rates[9] += (fwd_rates[381] - rev_rates[357]);
  //sp 5
  sp_rates[4] -= (fwd_rates[381] - rev_rates[357]);
  //sp 62
  sp_rates[61] -= (fwd_rates[381] - rev_rates[357]);
  //sp 63
  sp_rates[62] += (fwd_rates[381] - rev_rates[357]);

  //rxn 382
  //sp 62
  sp_rates[61] -= (fwd_rates[382] - rev_rates[358]);
  //sp 6
  sp_rates[5] -= (fwd_rates[382] - rev_rates[358]);
  //sp 63
  sp_rates[62] += (fwd_rates[382] - rev_rates[358]);
  //sp 8
  sp_rates[7] += (fwd_rates[382] - rev_rates[358]);

  //rxn 383
  //sp 63
  sp_rates[62] += (fwd_rates[383] - rev_rates[359]);
  //sp 62
  sp_rates[61] -= (fwd_rates[383] - rev_rates[359]);
  //sp 7
  sp_rates[6] += (fwd_rates[383] - rev_rates[359]);
  //sp 8
  sp_rates[7] -= (fwd_rates[383] - rev_rates[359]);

  //rxn 384
  //sp 9
  sp_rates[8] += (fwd_rates[384] - rev_rates[360]);
  //sp 10
  sp_rates[9] -= (fwd_rates[384] - rev_rates[360]);
  //sp 62
  sp_rates[61] -= (fwd_rates[384] - rev_rates[360]);
  //sp 63
  sp_rates[62] += (fwd_rates[384] - rev_rates[360]);

  //rxn 385
  //sp 63
  sp_rates[62] += (fwd_rates[385] - rev_rates[361]);
  //sp 14
  sp_rates[13] += (fwd_rates[385] - rev_rates[361]);
  //sp 62
  sp_rates[61] -= (fwd_rates[385] - rev_rates[361]);
  //sp 15
  sp_rates[14] -= (fwd_rates[385] - rev_rates[361]);

  //rxn 386
  //sp 63
  sp_rates[62] += (fwd_rates[386] - rev_rates[362]);
  //sp 22
  sp_rates[21] += (fwd_rates[386] - rev_rates[362]);
  //sp 62
  sp_rates[61] -= (fwd_rates[386] - rev_rates[362]);
  //sp 23
  sp_rates[22] -= (fwd_rates[386] - rev_rates[362]);

  //rxn 387
  //sp 20
  sp_rates[19] += (fwd_rates[387] - rev_rates[363]);
  //sp 21
  sp_rates[20] -= (fwd_rates[387] - rev_rates[363]);
  //sp 62
  sp_rates[61] -= (fwd_rates[387] - rev_rates[363]);
  //sp 63
  sp_rates[62] += (fwd_rates[387] - rev_rates[363]);

  //rxn 388
  //sp 63
  sp_rates[62] += (fwd_rates[388] - rev_rates[364]);
  //sp 30
  sp_rates[29] += (fwd_rates[388] - rev_rates[364]);
  //sp 62
  sp_rates[61] -= (fwd_rates[388] - rev_rates[364]);
  //sp 31
  sp_rates[30] -= (fwd_rates[388] - rev_rates[364]);

  //rxn 389
  //sp 33
  sp_rates[32] -= (fwd_rates[389] - rev_rates[365]);
  //sp 62
  sp_rates[61] -= (fwd_rates[389] - rev_rates[365]);
  //sp 63
  sp_rates[62] += (fwd_rates[389] - rev_rates[365]);
  //sp 32
  sp_rates[31] += (fwd_rates[389] - rev_rates[365]);

  //rxn 390
  //sp 65
  sp_rates[64] += (fwd_rates[390] - rev_rates[366]);
  //sp 3
  sp_rates[2] += (fwd_rates[390] - rev_rates[366]);
  //sp 4
  sp_rates[3] -= (fwd_rates[390] - rev_rates[366]);
  //sp 62
  sp_rates[61] -= (fwd_rates[390] - rev_rates[366]);

  //rxn 391
  //sp 65
  sp_rates[64] += (fwd_rates[391] - rev_rates[367]);
  //sp 62
  sp_rates[61] -= (fwd_rates[391] - rev_rates[367]);
  //sp 6
  sp_rates[5] -= (fwd_rates[391] - rev_rates[367]);
  //sp 8
  sp_rates[7] += (fwd_rates[391] - rev_rates[367]);

  //rxn 392
  //sp 65
  sp_rates[64] += (fwd_rates[392] - rev_rates[368]);
  //sp 62
  sp_rates[61] -= (fwd_rates[392] - rev_rates[368]);
  //sp 7
  sp_rates[6] += (fwd_rates[392] - rev_rates[368]);
  //sp 8
  sp_rates[7] -= (fwd_rates[392] - rev_rates[368]);

  //rxn 393
  //sp 9
  sp_rates[8] += (fwd_rates[393] - rev_rates[369]);
  //sp 10
  sp_rates[9] -= (fwd_rates[393] - rev_rates[369]);
  //sp 65
  sp_rates[64] += (fwd_rates[393] - rev_rates[369]);
  //sp 62
  sp_rates[61] -= (fwd_rates[393] - rev_rates[369]);

  //rxn 394
  //sp 65
  sp_rates[64] += (fwd_rates[394] - rev_rates[370]);
  //sp 10
  sp_rates[9] += (fwd_rates[394] - rev_rates[370]);
  //sp 5
  sp_rates[4] -= (fwd_rates[394] - rev_rates[370]);
  //sp 62
  sp_rates[61] -= (fwd_rates[394] - rev_rates[370]);

  //rxn 395
  //sp 65
  sp_rates[64] += (fwd_rates[395] - rev_rates[371]);
  //sp 14
  sp_rates[13] += (fwd_rates[395] - rev_rates[371]);
  //sp 62
  sp_rates[61] -= (fwd_rates[395] - rev_rates[371]);
  //sp 15
  sp_rates[14] -= (fwd_rates[395] - rev_rates[371]);

  //rxn 396
  //sp 3
  sp_rates[2] += (fwd_rates[396] - rev_rates[372]);
  //sp 4
  sp_rates[3] -= (fwd_rates[396] - rev_rates[372]);
  //sp 62
  sp_rates[61] -= (fwd_rates[396] - rev_rates[372]);
  //sp 64
  sp_rates[63] += (fwd_rates[396] - rev_rates[372]);

  //rxn 397
  //sp 3
  sp_rates[2] += (fwd_rates[397] - rev_rates[373]);
  //sp 4
  sp_rates[3] -= (fwd_rates[397] - rev_rates[373]);
  //sp 62
  sp_rates[61] -= (fwd_rates[397] - rev_rates[373]);
  //sp 64
  sp_rates[63] += (fwd_rates[397] - rev_rates[373]);

  //rxn 398
  //sp 10
  sp_rates[9] += (fwd_rates[398] - rev_rates[374]);
  //sp 5
  sp_rates[4] -= (fwd_rates[398] - rev_rates[374]);
  //sp 62
  sp_rates[61] -= (fwd_rates[398] - rev_rates[374]);
  //sp 64
  sp_rates[63] += (fwd_rates[398] - rev_rates[374]);

  //rxn 399
  //sp 64
  sp_rates[63] += (fwd_rates[399] - rev_rates[375]);
  //sp 62
  sp_rates[61] -= (fwd_rates[399] - rev_rates[375]);
  //sp 6
  sp_rates[5] -= (fwd_rates[399] - rev_rates[375]);
  //sp 8
  sp_rates[7] += (fwd_rates[399] - rev_rates[375]);

  //rxn 400
  //sp 64
  sp_rates[63] += (fwd_rates[400] - rev_rates[376]);
  //sp 62
  sp_rates[61] -= (fwd_rates[400] - rev_rates[376]);
  //sp 7
  sp_rates[6] += (fwd_rates[400] - rev_rates[376]);
  //sp 8
  sp_rates[7] -= (fwd_rates[400] - rev_rates[376]);

  //rxn 401
  //sp 9
  sp_rates[8] += (fwd_rates[401] - rev_rates[377]);
  //sp 10
  sp_rates[9] -= (fwd_rates[401] - rev_rates[377]);
  //sp 62
  sp_rates[61] -= (fwd_rates[401] - rev_rates[377]);
  //sp 64
  sp_rates[63] += (fwd_rates[401] - rev_rates[377]);

  //rxn 402
  //sp 14
  sp_rates[13] += (fwd_rates[402] - rev_rates[378]);
  //sp 62
  sp_rates[61] -= (fwd_rates[402] - rev_rates[378]);
  //sp 15
  sp_rates[14] -= (fwd_rates[402] - rev_rates[378]);
  //sp 64
  sp_rates[63] += (fwd_rates[402] - rev_rates[378]);

  //rxn 403
  //sp 26
  sp_rates[25] += (fwd_rates[403] - rev_rates[379]);
  //sp 62
  sp_rates[61] -= (fwd_rates[403] - rev_rates[379]);
  //sp 6
  sp_rates[5] -= (fwd_rates[403] - rev_rates[379]);
  //sp 31
  sp_rates[30] += (fwd_rates[403] - rev_rates[379]);

  //rxn 404
  //sp 4
  sp_rates[3] += fwd_rates[404];
  //sp 6
  sp_rates[5] -= fwd_rates[404];
  //sp 15
  sp_rates[14] += fwd_rates[404];
  //sp 48
  sp_rates[47] += fwd_rates[404];
  //sp 62
  sp_rates[61] -= fwd_rates[404];

  //rxn 405
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[405];
  //sp 62
  sp_rates[61] -= fwd_rates[405];
  //sp 70
  sp_rates[69] = fwd_rates[405];
  //sp 6
  sp_rates[5] -= fwd_rates[405];

  //rxn 406
  //sp 4
  sp_rates[3] -= (fwd_rates[406] - rev_rates[380]);
  //sp 61
  sp_rates[60] += (fwd_rates[406] - rev_rates[380]);
  //sp 62
  sp_rates[61] -= (fwd_rates[406] - rev_rates[380]);

  //rxn 407
  //sp 4
  sp_rates[3] -= (fwd_rates[407] - rev_rates[381]);
  //sp 61
  sp_rates[60] += (fwd_rates[407] - rev_rates[381]);
  //sp 62
  sp_rates[61] -= (fwd_rates[407] - rev_rates[381]);

  //rxn 408
  //sp 34
  sp_rates[33] += (fwd_rates[408] - rev_rates[382]);
  //sp 4
  sp_rates[3] -= (fwd_rates[408] - rev_rates[382]);
  //sp 62
  sp_rates[61] -= (fwd_rates[408] - rev_rates[382]);
  //sp 15
  sp_rates[14] += (fwd_rates[408] - rev_rates[382]);

  //rxn 409
  //sp 34
  sp_rates[33] += (fwd_rates[409] - rev_rates[383]);
  //sp 4
  sp_rates[3] -= (fwd_rates[409] - rev_rates[383]);
  //sp 62
  sp_rates[61] -= (fwd_rates[409] - rev_rates[383]);
  //sp 15
  sp_rates[14] += (fwd_rates[409] - rev_rates[383]);

  //rxn 410
  //sp 60
  sp_rates[59] += (fwd_rates[410] - rev_rates[384]);
  //sp 4
  sp_rates[3] -= (fwd_rates[410] - rev_rates[384]);
  //sp 62
  sp_rates[61] -= (fwd_rates[410] - rev_rates[384]);

  //rxn 411
  //sp 60
  sp_rates[59] += (fwd_rates[411] - rev_rates[385]);
  //sp 4
  sp_rates[3] -= (fwd_rates[411] - rev_rates[385]);
  //sp 62
  sp_rates[61] -= (fwd_rates[411] - rev_rates[385]);

  //rxn 412
  //sp 34
  sp_rates[33] -= (fwd_rates[412] - rev_rates[386]);
  //sp 61
  sp_rates[60] += (fwd_rates[412] - rev_rates[386]);
  //sp 15
  sp_rates[14] -= (fwd_rates[412] - rev_rates[386]);

  //rxn 413
  //sp 34
  sp_rates[33] -= (fwd_rates[413] - rev_rates[387]);
  //sp 61
  sp_rates[60] += (fwd_rates[413] - rev_rates[387]);
  //sp 15
  sp_rates[14] -= (fwd_rates[413] - rev_rates[387]);

  //rxn 414
  //sp 10
  sp_rates[9] -= (fwd_rates[414] - rev_rates[388]);
  //sp 60
  sp_rates[59] += (fwd_rates[414] - rev_rates[388]);
  //sp 5
  sp_rates[4] += (fwd_rates[414] - rev_rates[388]);
  //sp 62
  sp_rates[61] -= (fwd_rates[414] - rev_rates[388]);

  //rxn 415
  //sp 79
  sp_rates[78] = (fwd_rates[415] - rev_rates[389]);
  //sp 3
  sp_rates[2] += (fwd_rates[415] - rev_rates[389]);
  //sp 4
  sp_rates[3] -= (fwd_rates[415] - rev_rates[389]);
  //sp 63
  sp_rates[62] -= (fwd_rates[415] - rev_rates[389]);

  //rxn 416
  //sp 79
  sp_rates[78] += (fwd_rates[416] - rev_rates[390]);
  //sp 7
  sp_rates[6] += (fwd_rates[416] - rev_rates[390]);
  //sp 63
  sp_rates[62] -= (fwd_rates[416] - rev_rates[390]);
  //sp 8
  sp_rates[7] -= (fwd_rates[416] - rev_rates[390]);

  //rxn 417
  //sp 63
  sp_rates[62] -= (fwd_rates[417] - rev_rates[391]);
  //sp 79
  sp_rates[78] += (fwd_rates[417] - rev_rates[391]);
  //sp 14
  sp_rates[13] += (fwd_rates[417] - rev_rates[391]);
  //sp 15
  sp_rates[14] -= (fwd_rates[417] - rev_rates[391]);

  //rxn 418
  //sp 63
  sp_rates[62] -= (fwd_rates[418] - rev_rates[392]);
  //sp 79
  sp_rates[78] += (fwd_rates[418] - rev_rates[392]);
  //sp 30
  sp_rates[29] += (fwd_rates[418] - rev_rates[392]);
  //sp 31
  sp_rates[30] -= (fwd_rates[418] - rev_rates[392]);

  //rxn 419
  //sp 34
  sp_rates[33] += (fwd_rates[419] - rev_rates[393]);
  //sp 35
  sp_rates[34] -= (fwd_rates[419] - rev_rates[393]);
  //sp 63
  sp_rates[62] -= (fwd_rates[419] - rev_rates[393]);
  //sp 79
  sp_rates[78] += (fwd_rates[419] - rev_rates[393]);

  //rxn 420
  //sp 63
  sp_rates[62] += (fwd_rates[420] - rev_rates[394]);
  //sp 79
  sp_rates[78] -= 2.0 * (fwd_rates[420] - rev_rates[394]);
  //sp 80
  sp_rates[79] += (fwd_rates[420] - rev_rates[394]);

  //rxn 421
  //sp 3
  sp_rates[2] += (fwd_rates[421] - rev_rates[395]);
  //sp 4
  sp_rates[3] -= (fwd_rates[421] - rev_rates[395]);
  //sp 79
  sp_rates[78] += (fwd_rates[421] - rev_rates[395]);
  //sp 64
  sp_rates[63] -= (fwd_rates[421] - rev_rates[395]);

  //rxn 422
  //sp 79
  sp_rates[78] += (fwd_rates[422] - rev_rates[396]);
  //sp 14
  sp_rates[13] += (fwd_rates[422] - rev_rates[396]);
  //sp 15
  sp_rates[14] -= (fwd_rates[422] - rev_rates[396]);
  //sp 64
  sp_rates[63] -= (fwd_rates[422] - rev_rates[396]);

  //rxn 423
  //sp 3
  sp_rates[2] += (fwd_rates[423] - rev_rates[397]);
  //sp 4
  sp_rates[3] -= (fwd_rates[423] - rev_rates[397]);
  //sp 78
  sp_rates[77] = (fwd_rates[423] - rev_rates[397]);
  //sp 64
  sp_rates[63] -= (fwd_rates[423] - rev_rates[397]);

  //rxn 424
  //sp 78
  sp_rates[77] += (fwd_rates[424] - rev_rates[398]);
  //sp 14
  sp_rates[13] += (fwd_rates[424] - rev_rates[398]);
  //sp 15
  sp_rates[14] -= (fwd_rates[424] - rev_rates[398]);
  //sp 64
  sp_rates[63] -= (fwd_rates[424] - rev_rates[398]);

  //rxn 425
  //sp 65
  sp_rates[64] -= (fwd_rates[425] - rev_rates[399]);
  //sp 3
  sp_rates[2] += (fwd_rates[425] - rev_rates[399]);
  //sp 4
  sp_rates[3] -= (fwd_rates[425] - rev_rates[399]);
  //sp 78
  sp_rates[77] += (fwd_rates[425] - rev_rates[399]);

  //rxn 426
  //sp 65
  sp_rates[64] -= (fwd_rates[426] - rev_rates[400]);
  //sp 78
  sp_rates[77] += (fwd_rates[426] - rev_rates[400]);
  //sp 14
  sp_rates[13] += (fwd_rates[426] - rev_rates[400]);
  //sp 15
  sp_rates[14] -= (fwd_rates[426] - rev_rates[400]);

  //rxn 427
  //sp 79
  sp_rates[78] += (fwd_rates[427] - rev_rates[401]);
  //sp 62
  sp_rates[61] += (fwd_rates[427] - rev_rates[401]);
  //sp 63
  sp_rates[62] -= 2.0 * (fwd_rates[427] - rev_rates[401]);

  //rxn 428
  //sp 63
  sp_rates[62] -= (fwd_rates[428] - rev_rates[402]);
  //sp 34
  sp_rates[33] += (fwd_rates[428] - rev_rates[402]);
  //sp 62
  sp_rates[61] += (fwd_rates[428] - rev_rates[402]);
  //sp 31
  sp_rates[30] -= (fwd_rates[428] - rev_rates[402]);

  //rxn 429
  //sp 26
  sp_rates[25] -= (fwd_rates[429] - rev_rates[403]);
  //sp 11
  sp_rates[10] += (fwd_rates[429] - rev_rates[403]);
  //sp 62
  sp_rates[61] += (fwd_rates[429] - rev_rates[403]);
  //sp 63
  sp_rates[62] -= (fwd_rates[429] - rev_rates[403]);

  //rxn 430
  //sp 26
  sp_rates[25] -= (fwd_rates[430] - rev_rates[404]);
  //sp 11
  sp_rates[10] += (fwd_rates[430] - rev_rates[404]);
  //sp 62
  sp_rates[61] += (fwd_rates[430] - rev_rates[404]);
  //sp 64
  sp_rates[63] -= (fwd_rates[430] - rev_rates[404]);

  //rxn 431
  //sp 65
  sp_rates[64] -= (fwd_rates[431] - rev_rates[405]);
  //sp 26
  sp_rates[25] -= (fwd_rates[431] - rev_rates[405]);
  //sp 11
  sp_rates[10] += (fwd_rates[431] - rev_rates[405]);
  //sp 62
  sp_rates[61] += (fwd_rates[431] - rev_rates[405]);

  //rxn 432
  //sp 26
  sp_rates[25] += (fwd_rates[432] - rev_rates[406]);
  //sp 34
  sp_rates[33] += (fwd_rates[432] - rev_rates[406]);
  //sp 6
  sp_rates[5] -= (fwd_rates[432] - rev_rates[406]);
  //sp 64
  sp_rates[63] -= (fwd_rates[432] - rev_rates[406]);

  //rxn 433
  //sp 34
  sp_rates[33] += fwd_rates[433];
  //sp 4
  sp_rates[3] += fwd_rates[433];
  //sp 8
  sp_rates[7] -= fwd_rates[433];
  //sp 26
  sp_rates[25] += fwd_rates[433];
  //sp 64
  sp_rates[63] -= fwd_rates[433];

  //rxn 434
  //sp 34
  sp_rates[33] += fwd_rates[434];
  //sp 8
  sp_rates[7] += fwd_rates[434];
  //sp 10
  sp_rates[9] -= fwd_rates[434];
  //sp 26
  sp_rates[25] += fwd_rates[434];
  //sp 64
  sp_rates[63] -= fwd_rates[434];

  //rxn 435
  //sp 65
  sp_rates[64] -= (fwd_rates[435] - rev_rates[407]);
  //sp 6
  sp_rates[5] -= (fwd_rates[435] - rev_rates[407]);
  //sp 15
  sp_rates[14] += (fwd_rates[435] - rev_rates[407]);
  //sp 48
  sp_rates[47] += (fwd_rates[435] - rev_rates[407]);

  //rxn 436
  //sp 65
  sp_rates[64] -= fwd_rates[436];
  //sp 4
  sp_rates[3] += fwd_rates[436];
  //sp 8
  sp_rates[7] -= fwd_rates[436];
  //sp 15
  sp_rates[14] += fwd_rates[436];
  //sp 48
  sp_rates[47] += fwd_rates[436];

  //rxn 437
  //sp 65
  sp_rates[64] -= fwd_rates[437];
  //sp 8
  sp_rates[7] += fwd_rates[437];
  //sp 10
  sp_rates[9] -= fwd_rates[437];
  //sp 15
  sp_rates[14] += fwd_rates[437];
  //sp 48
  sp_rates[47] += fwd_rates[437];

  //rxn 438
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[438];
  //sp 85
  sp_rates[84] = fwd_rates[438];
  //sp 63
  sp_rates[62] -= fwd_rates[438];
  //sp 8
  sp_rates[7] -= fwd_rates[438];

  //rxn 439
  //sp 10
  sp_rates[9] += (fwd_rates[439] - rev_rates[408]);
  //sp 5
  sp_rates[4] -= (fwd_rates[439] - rev_rates[408]);
  //sp 63
  sp_rates[62] -= (fwd_rates[439] - rev_rates[408]);
  //sp 79
  sp_rates[78] += (fwd_rates[439] - rev_rates[408]);

  //rxn 440
  //sp 25
  sp_rates[24] += (fwd_rates[440] - rev_rates[409]);
  //sp 5
  sp_rates[4] -= (fwd_rates[440] - rev_rates[409]);
  //sp 46
  sp_rates[45] += (fwd_rates[440] - rev_rates[409]);
  //sp 63
  sp_rates[62] -= (fwd_rates[440] - rev_rates[409]);

  //rxn 441
  //sp 5
  sp_rates[4] -= (fwd_rates[441] - rev_rates[410]);
  //sp 85
  sp_rates[84] += (fwd_rates[441] - rev_rates[410]);
  //sp 63
  sp_rates[62] -= (fwd_rates[441] - rev_rates[410]);
  //sp 8
  sp_rates[7] += (fwd_rates[441] - rev_rates[410]);

  //rxn 442
  //sp 26
  sp_rates[25] += (fwd_rates[442] - rev_rates[411]);
  //sp 5
  sp_rates[4] -= (fwd_rates[442] - rev_rates[411]);
  //sp 45
  sp_rates[44] += (fwd_rates[442] - rev_rates[411]);
  //sp 64
  sp_rates[63] -= (fwd_rates[442] - rev_rates[411]);

  //rxn 443
  //sp 68
  sp_rates[67] = (fwd_rates[443] - rev_rates[412]);
  //sp 5
  sp_rates[4] -= (fwd_rates[443] - rev_rates[412]);
  //sp 6
  sp_rates[5] += (fwd_rates[443] - rev_rates[412]);
  //sp 64
  sp_rates[63] -= (fwd_rates[443] - rev_rates[412]);

  //rxn 444
  //sp 85
  sp_rates[84] += (fwd_rates[444] - rev_rates[413]);
  //sp 5
  sp_rates[4] -= (fwd_rates[444] - rev_rates[413]);
  //sp 8
  sp_rates[7] += (fwd_rates[444] - rev_rates[413]);
  //sp 64
  sp_rates[63] -= (fwd_rates[444] - rev_rates[413]);

  //rxn 445
  //sp 65
  sp_rates[64] -= (fwd_rates[445] - rev_rates[414]);
  //sp 5
  sp_rates[4] -= (fwd_rates[445] - rev_rates[414]);
  //sp 6
  sp_rates[5] += (fwd_rates[445] - rev_rates[414]);
  //sp 87
  sp_rates[86] = (fwd_rates[445] - rev_rates[414]);

  //rxn 446
  //sp 65
  sp_rates[64] -= (fwd_rates[446] - rev_rates[415]);
  //sp 25
  sp_rates[24] += (fwd_rates[446] - rev_rates[415]);
  //sp 5
  sp_rates[4] -= (fwd_rates[446] - rev_rates[415]);
  //sp 46
  sp_rates[45] += (fwd_rates[446] - rev_rates[415]);

  //rxn 447
  //sp 65
  sp_rates[64] -= (fwd_rates[447] - rev_rates[416]);
  //sp 10
  sp_rates[9] += (fwd_rates[447] - rev_rates[416]);
  //sp 5
  sp_rates[4] -= (fwd_rates[447] - rev_rates[416]);
  //sp 79
  sp_rates[78] += (fwd_rates[447] - rev_rates[416]);

  //rxn 448
  //sp 66
  sp_rates[65] = (fwd_rates[448] - rev_rates[417]);
  //sp 10
  sp_rates[9] -= (fwd_rates[448] - rev_rates[417]);
  //sp 63
  sp_rates[62] -= (fwd_rates[448] - rev_rates[417]);
  //sp 8
  sp_rates[7] += (fwd_rates[448] - rev_rates[417]);

  //rxn 449
  //sp 10
  sp_rates[9] -= (fwd_rates[449] - rev_rates[418]);
  //sp 63
  sp_rates[62] -= (fwd_rates[449] - rev_rates[418]);
  //sp 71
  sp_rates[70] = (fwd_rates[449] - rev_rates[418]);

  //rxn 450
  //sp 10
  sp_rates[9] -= (fwd_rates[450] - rev_rates[419]);
  //sp 85
  sp_rates[84] += (fwd_rates[450] - rev_rates[419]);
  //sp 63
  sp_rates[62] -= (fwd_rates[450] - rev_rates[419]);
  //sp 7
  sp_rates[6] += (fwd_rates[450] - rev_rates[419]);

  //rxn 451
  //sp 85
  sp_rates[84] += (fwd_rates[451] - rev_rates[420]);
  //sp 71
  sp_rates[70] -= (fwd_rates[451] - rev_rates[420]);
  //sp 7
  sp_rates[6] += (fwd_rates[451] - rev_rates[420]);

  //rxn 452
  //sp 66
  sp_rates[65] += (fwd_rates[452] - rev_rates[421]);
  //sp 71
  sp_rates[70] -= (fwd_rates[452] - rev_rates[421]);
  //sp 8
  sp_rates[7] += (fwd_rates[452] - rev_rates[421]);

  //rxn 453
  //sp 25
  sp_rates[24] += (fwd_rates[453] - rev_rates[422]);
  //sp 66
  sp_rates[65] -= (fwd_rates[453] - rev_rates[422]);
  //sp 35
  sp_rates[34] += (fwd_rates[453] - rev_rates[422]);

  //rxn 454
  //sp 66
  sp_rates[65] -= (fwd_rates[454] - rev_rates[423]);
  //sp 83
  sp_rates[82] = (fwd_rates[454] - rev_rates[423]);

  //rxn 455
  //sp 66
  sp_rates[65] -= (fwd_rates[455] - rev_rates[424]);
  //sp 4
  sp_rates[3] += (fwd_rates[455] - rev_rates[424]);
  //sp 85
  sp_rates[84] += (fwd_rates[455] - rev_rates[424]);

  //rxn 456
  //sp 66
  sp_rates[65] -= (fwd_rates[456] - rev_rates[425]);
  //sp 26
  sp_rates[25] += (fwd_rates[456] - rev_rates[425]);
  //sp 34
  sp_rates[33] += (fwd_rates[456] - rev_rates[425]);

  //rxn 457
  //sp 25
  sp_rates[24] += (fwd_rates[457] - rev_rates[426]);
  //sp 83
  sp_rates[82] -= (fwd_rates[457] - rev_rates[426]);
  //sp 35
  sp_rates[34] += (fwd_rates[457] - rev_rates[426]);

  //rxn 458
  //sp 83
  sp_rates[82] -= (fwd_rates[458] - rev_rates[427]);
  //sp 4
  sp_rates[3] += (fwd_rates[458] - rev_rates[427]);
  //sp 85
  sp_rates[84] += (fwd_rates[458] - rev_rates[427]);

  //rxn 459
  //sp 26
  sp_rates[25] += (fwd_rates[459] - rev_rates[428]);
  //sp 83
  sp_rates[82] -= (fwd_rates[459] - rev_rates[428]);
  //sp 34
  sp_rates[33] += (fwd_rates[459] - rev_rates[428]);

  //rxn 460
  //sp 25
  sp_rates[24] -= (fwd_rates[460] - rev_rates[429]);
  //sp 35
  sp_rates[34] -= (fwd_rates[460] - rev_rates[429]);
  //sp 4
  sp_rates[3] += (fwd_rates[460] - rev_rates[429]);
  //sp 85
  sp_rates[84] += (fwd_rates[460] - rev_rates[429]);

  //rxn 461
  //sp 25
  sp_rates[24] -= (fwd_rates[461] - rev_rates[430]);
  //sp 26
  sp_rates[25] += (fwd_rates[461] - rev_rates[430]);
  //sp 35
  sp_rates[34] -= (fwd_rates[461] - rev_rates[430]);
  //sp 34
  sp_rates[33] += (fwd_rates[461] - rev_rates[430]);

  //rxn 462
  //sp 66
  sp_rates[65] += (fwd_rates[462] - rev_rates[431]);
  //sp 21
  sp_rates[20] -= (fwd_rates[462] - rev_rates[431]);
  //sp 63
  sp_rates[62] -= (fwd_rates[462] - rev_rates[431]);
  //sp 23
  sp_rates[22] += (fwd_rates[462] - rev_rates[431]);

  //rxn 463
  //sp 21
  sp_rates[20] -= (fwd_rates[463] - rev_rates[432]);
  //sp 63
  sp_rates[62] -= (fwd_rates[463] - rev_rates[432]);
  //sp 69
  sp_rates[68] = (fwd_rates[463] - rev_rates[432]);

  //rxn 464
  //sp 66
  sp_rates[65] += (fwd_rates[464] - rev_rates[433]);
  //sp 69
  sp_rates[68] -= (fwd_rates[464] - rev_rates[433]);
  //sp 23
  sp_rates[22] += (fwd_rates[464] - rev_rates[433]);

  //rxn 465
  //sp 76
  sp_rates[75] = (fwd_rates[465] - rev_rates[434]);
  //sp 4
  sp_rates[3] += (fwd_rates[465] - rev_rates[434]);
  //sp 62
  sp_rates[61] -= (fwd_rates[465] - rev_rates[434]);
  //sp 8
  sp_rates[7] -= (fwd_rates[465] - rev_rates[434]);

  //rxn 466
  //sp 4
  sp_rates[3] += (fwd_rates[466] - rev_rates[435]);
  //sp 86
  sp_rates[85] += (fwd_rates[466] - rev_rates[435]);
  //sp 62
  sp_rates[61] -= (fwd_rates[466] - rev_rates[435]);
  //sp 8
  sp_rates[7] -= (fwd_rates[466] - rev_rates[435]);

  //rxn 467
  //sp 75
  sp_rates[74] = (fwd_rates[467] - rev_rates[436]);
  //sp 4
  sp_rates[3] += (fwd_rates[467] - rev_rates[436]);
  //sp 62
  sp_rates[61] -= (fwd_rates[467] - rev_rates[436]);
  //sp 8
  sp_rates[7] -= (fwd_rates[467] - rev_rates[436]);

  //rxn 468
  //sp 45
  sp_rates[44] += (fwd_rates[468] - rev_rates[437]);
  //sp 62
  sp_rates[61] -= (fwd_rates[468] - rev_rates[437]);
  //sp 15
  sp_rates[14] += (fwd_rates[468] - rev_rates[437]);
  //sp 8
  sp_rates[7] -= (fwd_rates[468] - rev_rates[437]);

  //rxn 469
  //sp 62
  sp_rates[61] -= (fwd_rates[469] - rev_rates[438]);
  //sp 72
  sp_rates[71] = (fwd_rates[469] - rev_rates[438]);
  //sp 8
  sp_rates[7] -= (fwd_rates[469] - rev_rates[438]);

  //rxn 470
  //sp 62
  sp_rates[61] -= (fwd_rates[470] - rev_rates[439]);
  //sp 72
  sp_rates[71] += (fwd_rates[470] - rev_rates[439]);
  //sp 8
  sp_rates[7] -= (fwd_rates[470] - rev_rates[439]);

  //rxn 471
  //sp 73
  sp_rates[72] = (fwd_rates[471] - rev_rates[440]);
  //sp 62
  sp_rates[61] -= (fwd_rates[471] - rev_rates[440]);
  //sp 8
  sp_rates[7] -= (fwd_rates[471] - rev_rates[440]);

  //rxn 472
  //sp 73
  sp_rates[72] += (fwd_rates[472] - rev_rates[441]);
  //sp 62
  sp_rates[61] -= (fwd_rates[472] - rev_rates[441]);
  //sp 8
  sp_rates[7] -= (fwd_rates[472] - rev_rates[441]);

  //rxn 473
  //sp 73
  sp_rates[72] += (fwd_rates[473] - rev_rates[442]);
  //sp 4
  sp_rates[3] -= (fwd_rates[473] - rev_rates[442]);
  //sp 86
  sp_rates[85] -= (fwd_rates[473] - rev_rates[442]);

  //rxn 474
  //sp 82
  sp_rates[81] = -(fwd_rates[474] - rev_rates[443]);
  //sp 4
  sp_rates[3] -= (fwd_rates[474] - rev_rates[443]);
  //sp 72
  sp_rates[71] += (fwd_rates[474] - rev_rates[443]);

  //rxn 475
  //sp 65
  sp_rates[64] -= (fwd_rates[475] - rev_rates[444]);
  //sp 86
  sp_rates[85] += (fwd_rates[475] - rev_rates[444]);
  //sp 8
  sp_rates[7] -= (fwd_rates[475] - rev_rates[444]);

  //rxn 476
  //sp 82
  sp_rates[81] += (fwd_rates[476] - rev_rates[445]);
  //sp 10
  sp_rates[9] += (fwd_rates[476] - rev_rates[445]);
  //sp 5
  sp_rates[4] -= (fwd_rates[476] - rev_rates[445]);
  //sp 72
  sp_rates[71] -= (fwd_rates[476] - rev_rates[445]);

  //rxn 477
  //sp 74
  sp_rates[73] = (fwd_rates[477] - rev_rates[446]);
  //sp 5
  sp_rates[4] -= (fwd_rates[477] - rev_rates[446]);
  //sp 72
  sp_rates[71] -= (fwd_rates[477] - rev_rates[446]);

  //rxn 478
  //sp 73
  sp_rates[72] -= (fwd_rates[478] - rev_rates[447]);
  //sp 74
  sp_rates[73] += (fwd_rates[478] - rev_rates[447]);
  //sp 5
  sp_rates[4] -= (fwd_rates[478] - rev_rates[447]);

  //rxn 479
  //sp 25
  sp_rates[24] += fwd_rates[479];
  //sp 74
  sp_rates[73] -= fwd_rates[479];
  //sp 45
  sp_rates[44] += fwd_rates[479];
  //sp 8
  sp_rates[7] += fwd_rates[479];

  //rxn 480
  //sp 82
  sp_rates[81] += (fwd_rates[480] - rev_rates[448]);
  //sp 75
  sp_rates[74] -= (fwd_rates[480] - rev_rates[448]);

  //rxn 481
  //sp 4
  sp_rates[3] += fwd_rates[481];
  //sp 5
  sp_rates[4] -= fwd_rates[481];
  //sp 10
  sp_rates[9] += fwd_rates[481];
  //sp 75
  sp_rates[74] -= fwd_rates[481];
  //sp 85
  sp_rates[84] += fwd_rates[481];

  //rxn 482
  //sp 4
  sp_rates[3] += fwd_rates[482];
  //sp 7
  sp_rates[6] += fwd_rates[482];
  //sp 8
  sp_rates[7] -= fwd_rates[482];
  //sp 75
  sp_rates[74] -= fwd_rates[482];
  //sp 85
  sp_rates[84] += fwd_rates[482];

  //rxn 483
  //sp 3
  sp_rates[2] += fwd_rates[483];
  //sp 75
  sp_rates[74] -= fwd_rates[483];
  //sp 85
  sp_rates[84] += fwd_rates[483];

  //rxn 484
  //sp 4
  sp_rates[3] += fwd_rates[484];
  //sp 6
  sp_rates[5] -= fwd_rates[484];
  //sp 8
  sp_rates[7] += fwd_rates[484];
  //sp 75
  sp_rates[74] -= fwd_rates[484];
  //sp 85
  sp_rates[84] += fwd_rates[484];

  //rxn 485
  //sp 4
  sp_rates[3] += fwd_rates[485];
  //sp 9
  sp_rates[8] += fwd_rates[485];
  //sp 10
  sp_rates[9] -= fwd_rates[485];
  //sp 75
  sp_rates[74] -= fwd_rates[485];
  //sp 85
  sp_rates[84] += fwd_rates[485];

  //rxn 486
  //sp 4
  sp_rates[3] += fwd_rates[486];
  //sp 75
  sp_rates[74] -= fwd_rates[486];
  //sp 14
  sp_rates[13] += fwd_rates[486];
  //sp 15
  sp_rates[14] -= fwd_rates[486];
  //sp 85
  sp_rates[84] += fwd_rates[486];

  //rxn 487
  //sp 4
  sp_rates[3] += fwd_rates[487];
  //sp 75
  sp_rates[74] -= fwd_rates[487];
  //sp 20
  sp_rates[19] += fwd_rates[487];
  //sp 21
  sp_rates[20] -= fwd_rates[487];
  //sp 85
  sp_rates[84] += fwd_rates[487];

  //rxn 488
  //sp 4
  sp_rates[3] += fwd_rates[488];
  //sp 75
  sp_rates[74] -= fwd_rates[488];
  //sp 85
  sp_rates[84] += fwd_rates[488];
  //sp 22
  sp_rates[21] += fwd_rates[488];
  //sp 23
  sp_rates[22] -= fwd_rates[488];

  //rxn 489
  //sp 82
  sp_rates[81] += (fwd_rates[489] - rev_rates[449]);
  //sp 75
  sp_rates[74] -= (fwd_rates[489] - rev_rates[449]);

  //rxn 490
  //sp 82
  sp_rates[81] += (fwd_rates[490] - rev_rates[450]);
  //sp 75
  sp_rates[74] -= (fwd_rates[490] - rev_rates[450]);

  //rxn 491
  //sp 78
  sp_rates[77] -= (fwd_rates[491] - rev_rates[451]);
  //sp 79
  sp_rates[78] += (fwd_rates[491] - rev_rates[451]);

  //rxn 492
  //sp 10
  sp_rates[9] += (fwd_rates[492] - rev_rates[452]);
  //sp 5
  sp_rates[4] -= (fwd_rates[492] - rev_rates[452]);
  //sp 78
  sp_rates[77] -= (fwd_rates[492] - rev_rates[452]);
  //sp 80
  sp_rates[79] += (fwd_rates[492] - rev_rates[452]);

  //rxn 493
  //sp 78
  sp_rates[77] -= (fwd_rates[493] - rev_rates[453]);
  //sp 6
  sp_rates[5] -= (fwd_rates[493] - rev_rates[453]);
  //sp 80
  sp_rates[79] += (fwd_rates[493] - rev_rates[453]);
  //sp 8
  sp_rates[7] += (fwd_rates[493] - rev_rates[453]);

  //rxn 494
  //sp 3
  sp_rates[2] += (fwd_rates[494] - rev_rates[454]);
  //sp 4
  sp_rates[3] -= (fwd_rates[494] - rev_rates[454]);
  //sp 78
  sp_rates[77] -= (fwd_rates[494] - rev_rates[454]);
  //sp 80
  sp_rates[79] += (fwd_rates[494] - rev_rates[454]);

  //rxn 495
  //sp 80
  sp_rates[79] += (fwd_rates[495] - rev_rates[455]);
  //sp 78
  sp_rates[77] -= (fwd_rates[495] - rev_rates[455]);
  //sp 7
  sp_rates[6] += (fwd_rates[495] - rev_rates[455]);
  //sp 8
  sp_rates[7] -= (fwd_rates[495] - rev_rates[455]);

  //rxn 496
  //sp 8
  sp_rates[7] -= fwd_rates[496];
  //sp 11
  sp_rates[10] += 0.15 * fwd_rates[496];
  //sp 78
  sp_rates[77] -= fwd_rates[496];
  //sp 87
  sp_rates[86] += 0.85 * fwd_rates[496];
  //sp 31
  sp_rates[30] += 0.15 * fwd_rates[496];

  //rxn 497
  //sp 9
  sp_rates[8] += (fwd_rates[497] - rev_rates[456]);
  //sp 10
  sp_rates[9] -= (fwd_rates[497] - rev_rates[456]);
  //sp 78
  sp_rates[77] -= (fwd_rates[497] - rev_rates[456]);
  //sp 80
  sp_rates[79] += (fwd_rates[497] - rev_rates[456]);

  //rxn 498
  //sp 14
  sp_rates[13] += (fwd_rates[498] - rev_rates[457]);
  //sp 78
  sp_rates[77] -= (fwd_rates[498] - rev_rates[457]);
  //sp 15
  sp_rates[14] -= (fwd_rates[498] - rev_rates[457]);
  //sp 80
  sp_rates[79] += (fwd_rates[498] - rev_rates[457]);

  //rxn 499
  //sp 20
  sp_rates[19] += (fwd_rates[499] - rev_rates[458]);
  //sp 21
  sp_rates[20] -= (fwd_rates[499] - rev_rates[458]);
  //sp 78
  sp_rates[77] -= (fwd_rates[499] - rev_rates[458]);
  //sp 80
  sp_rates[79] += (fwd_rates[499] - rev_rates[458]);

  //rxn 500
  //sp 36
  sp_rates[35] += (fwd_rates[500] - rev_rates[459]);
  //sp 37
  sp_rates[36] -= (fwd_rates[500] - rev_rates[459]);
  //sp 78
  sp_rates[77] -= (fwd_rates[500] - rev_rates[459]);
  //sp 80
  sp_rates[79] += (fwd_rates[500] - rev_rates[459]);

  //rxn 501
  //sp 34
  sp_rates[33] += (fwd_rates[501] - rev_rates[460]);
  //sp 35
  sp_rates[34] -= (fwd_rates[501] - rev_rates[460]);
  //sp 78
  sp_rates[77] -= (fwd_rates[501] - rev_rates[460]);
  //sp 80
  sp_rates[79] += (fwd_rates[501] - rev_rates[460]);

  //rxn 502
  //sp 62
  sp_rates[61] += (fwd_rates[502] - rev_rates[461]);
  //sp 78
  sp_rates[77] -= (fwd_rates[502] - rev_rates[461]);
  //sp 63
  sp_rates[62] -= (fwd_rates[502] - rev_rates[461]);
  //sp 80
  sp_rates[79] += (fwd_rates[502] - rev_rates[461]);

  //rxn 503
  //sp 3
  sp_rates[2] += (fwd_rates[503] - rev_rates[462]);
  //sp 4
  sp_rates[3] -= (fwd_rates[503] - rev_rates[462]);
  //sp 79
  sp_rates[78] -= (fwd_rates[503] - rev_rates[462]);
  //sp 80
  sp_rates[79] += (fwd_rates[503] - rev_rates[462]);

  //rxn 504
  //sp 10
  sp_rates[9] += (fwd_rates[504] - rev_rates[463]);
  //sp 5
  sp_rates[4] -= (fwd_rates[504] - rev_rates[463]);
  //sp 79
  sp_rates[78] -= (fwd_rates[504] - rev_rates[463]);
  //sp 80
  sp_rates[79] += (fwd_rates[504] - rev_rates[463]);

  //rxn 505
  //sp 7
  sp_rates[6] += (fwd_rates[505] - rev_rates[464]);
  //sp 80
  sp_rates[79] += (fwd_rates[505] - rev_rates[464]);
  //sp 79
  sp_rates[78] -= (fwd_rates[505] - rev_rates[464]);
  //sp 8
  sp_rates[7] -= (fwd_rates[505] - rev_rates[464]);

  //rxn 506
  //sp 79
  sp_rates[78] -= (fwd_rates[506] - rev_rates[465]);
  //sp 14
  sp_rates[13] += (fwd_rates[506] - rev_rates[465]);
  //sp 15
  sp_rates[14] -= (fwd_rates[506] - rev_rates[465]);
  //sp 80
  sp_rates[79] += (fwd_rates[506] - rev_rates[465]);

  //rxn 507
  //sp 9
  sp_rates[8] += (fwd_rates[507] - rev_rates[466]);
  //sp 10
  sp_rates[9] -= (fwd_rates[507] - rev_rates[466]);
  //sp 79
  sp_rates[78] -= (fwd_rates[507] - rev_rates[466]);
  //sp 80
  sp_rates[79] += (fwd_rates[507] - rev_rates[466]);

  //rxn 508
  //sp 20
  sp_rates[19] += (fwd_rates[508] - rev_rates[467]);
  //sp 21
  sp_rates[20] -= (fwd_rates[508] - rev_rates[467]);
  //sp 79
  sp_rates[78] -= (fwd_rates[508] - rev_rates[467]);
  //sp 80
  sp_rates[79] += (fwd_rates[508] - rev_rates[467]);

  //rxn 509
  //sp 79
  sp_rates[78] -= (fwd_rates[509] - rev_rates[468]);
  //sp 62
  sp_rates[61] += (fwd_rates[509] - rev_rates[468]);
  //sp 63
  sp_rates[62] -= (fwd_rates[509] - rev_rates[468]);
  //sp 80
  sp_rates[79] += (fwd_rates[509] - rev_rates[468]);

  //rxn 510
  //sp 78
  sp_rates[77] += (fwd_rates[510] - rev_rates[469]);
  //sp 79
  sp_rates[78] -= (fwd_rates[510] - rev_rates[469]);

  //rxn 511
  //sp 78
  sp_rates[77] += (fwd_rates[511] - rev_rates[470]);
  //sp 79
  sp_rates[78] -= (fwd_rates[511] - rev_rates[470]);

  //rxn 512
  //sp 63
  sp_rates[62] += (fwd_rates[512] - rev_rates[471]);
  //sp 4
  sp_rates[3] -= (fwd_rates[512] - rev_rates[471]);
  //sp 79
  sp_rates[78] -= (fwd_rates[512] - rev_rates[471]);

  //rxn 513
  //sp 63
  sp_rates[62] += (fwd_rates[513] - rev_rates[472]);
  //sp 4
  sp_rates[3] -= (fwd_rates[513] - rev_rates[472]);
  //sp 79
  sp_rates[78] -= (fwd_rates[513] - rev_rates[472]);

  //rxn 514
  //sp 4
  sp_rates[3] -= (fwd_rates[514] - rev_rates[473]);
  //sp 79
  sp_rates[78] -= (fwd_rates[514] - rev_rates[473]);
  //sp 64
  sp_rates[63] += (fwd_rates[514] - rev_rates[473]);

  //rxn 515
  //sp 4
  sp_rates[3] -= (fwd_rates[515] - rev_rates[474]);
  //sp 79
  sp_rates[78] -= (fwd_rates[515] - rev_rates[474]);
  //sp 64
  sp_rates[63] += (fwd_rates[515] - rev_rates[474]);

  //rxn 516
  //sp 65
  sp_rates[64] += (fwd_rates[516] - rev_rates[475]);
  //sp 4
  sp_rates[3] -= (fwd_rates[516] - rev_rates[475]);
  //sp 79
  sp_rates[78] -= (fwd_rates[516] - rev_rates[475]);

  //rxn 517
  //sp 65
  sp_rates[64] += (fwd_rates[517] - rev_rates[476]);
  //sp 4
  sp_rates[3] -= (fwd_rates[517] - rev_rates[476]);
  //sp 79
  sp_rates[78] -= (fwd_rates[517] - rev_rates[476]);

  //rxn 518
  //sp 36
  sp_rates[35] += (fwd_rates[518] - rev_rates[477]);
  //sp 4
  sp_rates[3] -= (fwd_rates[518] - rev_rates[477]);
  //sp 79
  sp_rates[78] -= (fwd_rates[518] - rev_rates[477]);
  //sp 15
  sp_rates[14] += (fwd_rates[518] - rev_rates[477]);

  //rxn 519
  //sp 36
  sp_rates[35] += (fwd_rates[519] - rev_rates[478]);
  //sp 4
  sp_rates[3] -= (fwd_rates[519] - rev_rates[478]);
  //sp 79
  sp_rates[78] -= (fwd_rates[519] - rev_rates[478]);
  //sp 15
  sp_rates[14] += (fwd_rates[519] - rev_rates[478]);

  //rxn 520
  //sp 65
  sp_rates[64] += (fwd_rates[520] - rev_rates[479]);
  //sp 4
  sp_rates[3] -= (fwd_rates[520] - rev_rates[479]);
  //sp 78
  sp_rates[77] -= (fwd_rates[520] - rev_rates[479]);

  //rxn 521
  //sp 65
  sp_rates[64] += (fwd_rates[521] - rev_rates[480]);
  //sp 4
  sp_rates[3] -= (fwd_rates[521] - rev_rates[480]);
  //sp 78
  sp_rates[77] -= (fwd_rates[521] - rev_rates[480]);

  //rxn 522
  //sp 4
  sp_rates[3] -= (fwd_rates[522] - rev_rates[481]);
  //sp 78
  sp_rates[77] -= (fwd_rates[522] - rev_rates[481]);
  //sp 64
  sp_rates[63] += (fwd_rates[522] - rev_rates[481]);

  //rxn 523
  //sp 4
  sp_rates[3] -= (fwd_rates[523] - rev_rates[482]);
  //sp 78
  sp_rates[77] -= (fwd_rates[523] - rev_rates[482]);
  //sp 64
  sp_rates[63] += (fwd_rates[523] - rev_rates[482]);

  //rxn 524
  //sp 36
  sp_rates[35] += (fwd_rates[524] - rev_rates[483]);
  //sp 4
  sp_rates[3] -= (fwd_rates[524] - rev_rates[483]);
  //sp 78
  sp_rates[77] -= (fwd_rates[524] - rev_rates[483]);
  //sp 15
  sp_rates[14] += (fwd_rates[524] - rev_rates[483]);

  //rxn 525
  //sp 4
  sp_rates[3] -= (fwd_rates[525] - rev_rates[484]);
  //sp 78
  sp_rates[77] -= (fwd_rates[525] - rev_rates[484]);
  //sp 63
  sp_rates[62] += (fwd_rates[525] - rev_rates[484]);

  //rxn 526
  //sp 65
  sp_rates[64] += (fwd_rates[526] - rev_rates[485]);
  //sp 63
  sp_rates[62] -= (fwd_rates[526] - rev_rates[485]);

  //rxn 527
  //sp 63
  sp_rates[62] -= (fwd_rates[527] - rev_rates[486]);
  //sp 64
  sp_rates[63] += (fwd_rates[527] - rev_rates[486]);

  //rxn 528
  //sp 65
  sp_rates[64] += (fwd_rates[528] - rev_rates[487]);
  //sp 36
  sp_rates[35] -= (fwd_rates[528] - rev_rates[487]);
  //sp 15
  sp_rates[14] -= (fwd_rates[528] - rev_rates[487]);

  //rxn 529
  //sp 65
  sp_rates[64] -= (fwd_rates[529] - rev_rates[488]);
  //sp 64
  sp_rates[63] += (fwd_rates[529] - rev_rates[488]);

  //rxn 530
  //sp 63
  sp_rates[62] += (fwd_rates[530] - rev_rates[489]);
  //sp 36
  sp_rates[35] -= (fwd_rates[530] - rev_rates[489]);
  //sp 15
  sp_rates[14] -= (fwd_rates[530] - rev_rates[489]);

  //rxn 531
  //sp 36
  sp_rates[35] -= (fwd_rates[531] - rev_rates[490]);
  //sp 15
  sp_rates[14] -= (fwd_rates[531] - rev_rates[490]);
  //sp 64
  sp_rates[63] += (fwd_rates[531] - rev_rates[490]);

  //rxn 532
  //sp 36
  sp_rates[35] -= (fwd_rates[532] - rev_rates[491]);
  //sp 15
  sp_rates[14] -= (fwd_rates[532] - rev_rates[491]);
  //sp 64
  sp_rates[63] += (fwd_rates[532] - rev_rates[491]);

  //rxn 533
  //sp 49
  sp_rates[48] += (fwd_rates[533] - rev_rates[492]);
  //sp 78
  sp_rates[77] -= (fwd_rates[533] - rev_rates[492]);
  //sp 6
  sp_rates[5] -= (fwd_rates[533] - rev_rates[492]);
  //sp 15
  sp_rates[14] += (fwd_rates[533] - rev_rates[492]);

  //rxn 534
  //sp 34
  sp_rates[33] += (fwd_rates[534] - rev_rates[493]);
  //sp 11
  sp_rates[10] += (fwd_rates[534] - rev_rates[493]);
  //sp 78
  sp_rates[77] -= (fwd_rates[534] - rev_rates[493]);
  //sp 6
  sp_rates[5] -= (fwd_rates[534] - rev_rates[493]);

  //rxn 535
  //sp 26
  sp_rates[25] += (fwd_rates[535] - rev_rates[494]);
  //sp 35
  sp_rates[34] += (fwd_rates[535] - rev_rates[494]);
  //sp 78
  sp_rates[77] -= (fwd_rates[535] - rev_rates[494]);
  //sp 6
  sp_rates[5] -= (fwd_rates[535] - rev_rates[494]);

  //rxn 536
  //sp 34
  sp_rates[33] += (fwd_rates[536] - rev_rates[495]);
  //sp 11
  sp_rates[10] += (fwd_rates[536] - rev_rates[495]);
  //sp 6
  sp_rates[5] -= (fwd_rates[536] - rev_rates[495]);
  //sp 79
  sp_rates[78] -= (fwd_rates[536] - rev_rates[495]);

  //rxn 537
  //sp 25
  sp_rates[24] += (fwd_rates[537] - rev_rates[496]);
  //sp 36
  sp_rates[35] += (fwd_rates[537] - rev_rates[496]);
  //sp 6
  sp_rates[5] -= (fwd_rates[537] - rev_rates[496]);
  //sp 79
  sp_rates[78] -= (fwd_rates[537] - rev_rates[496]);

  //rxn 538
  //sp 34
  sp_rates[33] += fwd_rates[538];
  //sp 8
  sp_rates[7] += fwd_rates[538];
  //sp 10
  sp_rates[9] -= fwd_rates[538];
  //sp 11
  sp_rates[10] += fwd_rates[538];
  //sp 78
  sp_rates[77] -= fwd_rates[538];

  //rxn 539
  //sp 34
  sp_rates[33] += fwd_rates[539];
  //sp 8
  sp_rates[7] += fwd_rates[539];
  //sp 10
  sp_rates[9] -= fwd_rates[539];
  //sp 11
  sp_rates[10] += fwd_rates[539];
  //sp 79
  sp_rates[78] -= fwd_rates[539];

  //rxn 540
  //sp 8
  sp_rates[7] += fwd_rates[540];
  //sp 10
  sp_rates[9] -= fwd_rates[540];
  //sp 79
  sp_rates[78] -= fwd_rates[540];
  //sp 16
  sp_rates[15] += fwd_rates[540];
  //sp 48
  sp_rates[47] += fwd_rates[540];

  //rxn 541
  //sp 36
  sp_rates[35] += (fwd_rates[541] - rev_rates[497]);
  //sp 37
  sp_rates[36] -= (fwd_rates[541] - rev_rates[497]);
  //sp 79
  sp_rates[78] -= (fwd_rates[541] - rev_rates[497]);
  //sp 80
  sp_rates[79] += (fwd_rates[541] - rev_rates[497]);

  //rxn 542
  //sp 77
  sp_rates[76] = (fwd_rates[542] - rev_rates[498]);
  //sp 79
  sp_rates[78] -= (fwd_rates[542] - rev_rates[498]);
  //sp 8
  sp_rates[7] -= (fwd_rates[542] - rev_rates[498]);

  //rxn 543
  //sp 36
  sp_rates[35] += (fwd_rates[543] - rev_rates[499]);
  //sp 77
  sp_rates[76] -= (fwd_rates[543] - rev_rates[499]);
  //sp 24
  sp_rates[23] += (fwd_rates[543] - rev_rates[499]);

  //rxn 544
  //sp 25
  sp_rates[24] += (fwd_rates[544] - rev_rates[500]);
  //sp 35
  sp_rates[34] += (fwd_rates[544] - rev_rates[500]);
  //sp 77
  sp_rates[76] -= (fwd_rates[544] - rev_rates[500]);

  //rxn 545
  //sp 4
  sp_rates[3] += (fwd_rates[545] - rev_rates[501]);
  //sp 77
  sp_rates[76] -= (fwd_rates[545] - rev_rates[501]);
  //sp 85
  sp_rates[84] += (fwd_rates[545] - rev_rates[501]);

  //rxn 546
  //sp 26
  sp_rates[25] += (fwd_rates[546] - rev_rates[502]);
  //sp 77
  sp_rates[76] -= (fwd_rates[546] - rev_rates[502]);
  //sp 34
  sp_rates[33] += (fwd_rates[546] - rev_rates[502]);

  //rxn 547
  //sp 4
  sp_rates[3] += (fwd_rates[547] - rev_rates[503]);
  //sp 68
  sp_rates[67] -= (fwd_rates[547] - rev_rates[503]);
  //sp 85
  sp_rates[84] += (fwd_rates[547] - rev_rates[503]);

  //rxn 548
  //sp 4
  sp_rates[3] += (fwd_rates[548] - rev_rates[504]);
  //sp 68
  sp_rates[67] -= (fwd_rates[548] - rev_rates[504]);
  //sp 70
  sp_rates[69] += (fwd_rates[548] - rev_rates[504]);

  //rxn 549
  //sp 25
  sp_rates[24] += (fwd_rates[549] - rev_rates[505]);
  //sp 37
  sp_rates[36] += (fwd_rates[549] - rev_rates[505]);
  //sp 6
  sp_rates[5] -= (fwd_rates[549] - rev_rates[505]);
  //sp 80
  sp_rates[79] -= (fwd_rates[549] - rev_rates[505]);

  //rxn 550
  //sp 35
  sp_rates[34] += fwd_rates[550];
  //sp 8
  sp_rates[7] += fwd_rates[550];
  //sp 10
  sp_rates[9] -= fwd_rates[550];
  //sp 11
  sp_rates[10] += fwd_rates[550];
  //sp 80
  sp_rates[79] -= fwd_rates[550];

  //rxn 551
  //sp 26
  sp_rates[25] -= (fwd_rates[551] - rev_rates[506]);
  //sp 11
  sp_rates[10] += (fwd_rates[551] - rev_rates[506]);
  //sp 79
  sp_rates[78] += (fwd_rates[551] - rev_rates[506]);
  //sp 80
  sp_rates[79] -= (fwd_rates[551] - rev_rates[506]);

  //rxn 552
  //sp 26
  sp_rates[25] -= (fwd_rates[552] - rev_rates[507]);
  //sp 11
  sp_rates[10] += (fwd_rates[552] - rev_rates[507]);
  //sp 78
  sp_rates[77] += (fwd_rates[552] - rev_rates[507]);
  //sp 80
  sp_rates[79] -= (fwd_rates[552] - rev_rates[507]);

  //rxn 553
  //sp 37
  sp_rates[36] -= (fwd_rates[553] - rev_rates[508]);
  //sp 80
  sp_rates[79] += (fwd_rates[553] - rev_rates[508]);
  //sp 31
  sp_rates[30] -= (fwd_rates[553] - rev_rates[508]);
  //sp 15
  sp_rates[14] += (fwd_rates[553] - rev_rates[508]);

  //rxn 554
  //sp 26
  sp_rates[25] += (fwd_rates[554] - rev_rates[509]);
  //sp 5
  sp_rates[4] -= (fwd_rates[554] - rev_rates[509]);
  //sp 48
  sp_rates[47] += (fwd_rates[554] - rev_rates[509]);
  //sp 80
  sp_rates[79] -= (fwd_rates[554] - rev_rates[509]);

  //rxn 555
  //sp 81
  sp_rates[80] = (fwd_rates[555] - rev_rates[510]);
  //sp 3
  sp_rates[2] += (fwd_rates[555] - rev_rates[510]);
  //sp 4
  sp_rates[3] -= (fwd_rates[555] - rev_rates[510]);
  //sp 80
  sp_rates[79] -= (fwd_rates[555] - rev_rates[510]);

  //rxn 556
  //sp 81
  sp_rates[80] += (fwd_rates[556] - rev_rates[511]);
  //sp 7
  sp_rates[6] += (fwd_rates[556] - rev_rates[511]);
  //sp 80
  sp_rates[79] -= (fwd_rates[556] - rev_rates[511]);
  //sp 8
  sp_rates[7] -= (fwd_rates[556] - rev_rates[511]);

  //rxn 557
  //sp 25
  sp_rates[24] += (fwd_rates[557] - rev_rates[512]);
  //sp 36
  sp_rates[35] += (fwd_rates[557] - rev_rates[512]);
  //sp 80
  sp_rates[79] -= (fwd_rates[557] - rev_rates[512]);
  //sp 8
  sp_rates[7] -= (fwd_rates[557] - rev_rates[512]);

  //rxn 558
  //sp 26
  sp_rates[25] += (fwd_rates[558] - rev_rates[513]);
  //sp 35
  sp_rates[34] += (fwd_rates[558] - rev_rates[513]);
  //sp 80
  sp_rates[79] -= (fwd_rates[558] - rev_rates[513]);
  //sp 8
  sp_rates[7] -= (fwd_rates[558] - rev_rates[513]);

  //rxn 559
  //sp 34
  sp_rates[33] += (fwd_rates[559] - rev_rates[514]);
  //sp 11
  sp_rates[10] += (fwd_rates[559] - rev_rates[514]);
  //sp 80
  sp_rates[79] -= (fwd_rates[559] - rev_rates[514]);
  //sp 8
  sp_rates[7] -= (fwd_rates[559] - rev_rates[514]);

  //rxn 560
  //sp 4
  sp_rates[3] += (fwd_rates[560] - rev_rates[515]);
  //sp 109
  sp_rates[108] = (fwd_rates[560] - rev_rates[515]);
  //sp 80
  sp_rates[79] -= 2.0 * (fwd_rates[560] - rev_rates[515]);

  //rxn 561
  //sp 107
  sp_rates[106] = (fwd_rates[561] - rev_rates[516]);
  //sp 80
  sp_rates[79] -= 2.0 * (fwd_rates[561] - rev_rates[516]);

  //rxn 562
  //sp 108
  sp_rates[107] = (fwd_rates[562] - rev_rates[517]);
  //sp 80
  sp_rates[79] -= 2.0 * (fwd_rates[562] - rev_rates[517]);

  //rxn 563
  //sp 108
  sp_rates[107] += fwd_rates[563];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[563];
  //sp 63
  sp_rates[62] -= fwd_rates[563];
  //sp 80
  sp_rates[79] -= fwd_rates[563];

  //rxn 564
  //sp 107
  sp_rates[106] += (fwd_rates[564] - rev_rates[518]);
  //sp 4
  sp_rates[3] += (fwd_rates[564] - rev_rates[518]);
  //sp 79
  sp_rates[78] -= (fwd_rates[564] - rev_rates[518]);
  //sp 80
  sp_rates[79] -= (fwd_rates[564] - rev_rates[518]);

  //rxn 565
  //sp 12
  sp_rates[11] += (fwd_rates[565] - rev_rates[519]);
  //sp 70
  sp_rates[69] -= (fwd_rates[565] - rev_rates[519]);
  //sp 31
  sp_rates[30] += (fwd_rates[565] - rev_rates[519]);
  //sp 8
  sp_rates[7] -= (fwd_rates[565] - rev_rates[519]);

  //rxn 566
  //sp 41
  sp_rates[40] += (fwd_rates[566] - rev_rates[520]);
  //sp 11
  sp_rates[10] += (fwd_rates[566] - rev_rates[520]);
  //sp 70
  sp_rates[69] -= (fwd_rates[566] - rev_rates[520]);
  //sp 8
  sp_rates[7] -= (fwd_rates[566] - rev_rates[520]);

  //rxn 567
  //sp 11
  sp_rates[10] += (fwd_rates[567] - rev_rates[521]);
  //sp 4
  sp_rates[3] -= (fwd_rates[567] - rev_rates[521]);
  //sp 70
  sp_rates[69] -= (fwd_rates[567] - rev_rates[521]);
  //sp 31
  sp_rates[30] += (fwd_rates[567] - rev_rates[521]);

  //rxn 568
  //sp 11
  sp_rates[10] += (fwd_rates[568] - rev_rates[522]);
  //sp 70
  sp_rates[69] -= (fwd_rates[568] - rev_rates[522]);
  //sp 6
  sp_rates[5] -= (fwd_rates[568] - rev_rates[522]);
  //sp 45
  sp_rates[44] += (fwd_rates[568] - rev_rates[522]);

  //rxn 569
  //sp 46
  sp_rates[45] += (fwd_rates[569] - rev_rates[523]);
  //sp 86
  sp_rates[85] -= (fwd_rates[569] - rev_rates[523]);
  //sp 15
  sp_rates[14] += (fwd_rates[569] - rev_rates[523]);

  //rxn 570
  //sp 4
  sp_rates[3] -= (fwd_rates[570] - rev_rates[524]);
  //sp 86
  sp_rates[85] += (fwd_rates[570] - rev_rates[524]);
  //sp 87
  sp_rates[86] -= (fwd_rates[570] - rev_rates[524]);

  //rxn 571
  //sp 87
  sp_rates[86] += (fwd_rates[571] - rev_rates[525]);
  //sp 86
  sp_rates[85] -= (fwd_rates[571] - rev_rates[525]);
  //sp 7
  sp_rates[6] += (fwd_rates[571] - rev_rates[525]);
  //sp 8
  sp_rates[7] -= (fwd_rates[571] - rev_rates[525]);

  //rxn 572
  //sp 3
  sp_rates[2] += (fwd_rates[572] - rev_rates[526]);
  //sp 4
  sp_rates[3] -= (fwd_rates[572] - rev_rates[526]);
  //sp 86
  sp_rates[85] -= (fwd_rates[572] - rev_rates[526]);
  //sp 87
  sp_rates[86] += (fwd_rates[572] - rev_rates[526]);

  //rxn 573
  //sp 86
  sp_rates[85] -= (fwd_rates[573] - rev_rates[527]);
  //sp 6
  sp_rates[5] -= (fwd_rates[573] - rev_rates[527]);
  //sp 87
  sp_rates[86] += (fwd_rates[573] - rev_rates[527]);
  //sp 8
  sp_rates[7] += (fwd_rates[573] - rev_rates[527]);

  //rxn 574
  //sp 87
  sp_rates[86] += (fwd_rates[574] - rev_rates[528]);
  //sp 14
  sp_rates[13] += (fwd_rates[574] - rev_rates[528]);
  //sp 86
  sp_rates[85] -= (fwd_rates[574] - rev_rates[528]);
  //sp 15
  sp_rates[14] -= (fwd_rates[574] - rev_rates[528]);

  //rxn 575
  //sp 87
  sp_rates[86] += (fwd_rates[575] - rev_rates[529]);
  //sp 22
  sp_rates[21] += (fwd_rates[575] - rev_rates[529]);
  //sp 86
  sp_rates[85] -= (fwd_rates[575] - rev_rates[529]);
  //sp 23
  sp_rates[22] -= (fwd_rates[575] - rev_rates[529]);

  //rxn 576
  //sp 10
  sp_rates[9] += (fwd_rates[576] - rev_rates[530]);
  //sp 5
  sp_rates[4] -= (fwd_rates[576] - rev_rates[530]);
  //sp 86
  sp_rates[85] -= (fwd_rates[576] - rev_rates[530]);
  //sp 87
  sp_rates[86] += (fwd_rates[576] - rev_rates[530]);

  //rxn 577
  //sp 9
  sp_rates[8] += (fwd_rates[577] - rev_rates[531]);
  //sp 10
  sp_rates[9] -= (fwd_rates[577] - rev_rates[531]);
  //sp 86
  sp_rates[85] -= (fwd_rates[577] - rev_rates[531]);
  //sp 87
  sp_rates[86] += (fwd_rates[577] - rev_rates[531]);

  //rxn 578
  //sp 20
  sp_rates[19] += (fwd_rates[578] - rev_rates[532]);
  //sp 21
  sp_rates[20] -= (fwd_rates[578] - rev_rates[532]);
  //sp 86
  sp_rates[85] -= (fwd_rates[578] - rev_rates[532]);
  //sp 87
  sp_rates[86] += (fwd_rates[578] - rev_rates[532]);

  //rxn 579
  //sp 26
  sp_rates[25] -= (fwd_rates[579] - rev_rates[533]);
  //sp 35
  sp_rates[34] -= (fwd_rates[579] - rev_rates[533]);
  //sp 85
  sp_rates[84] += (fwd_rates[579] - rev_rates[533]);

  //rxn 580
  //sp 31
  sp_rates[30] += 2.0 * (fwd_rates[580] - rev_rates[534]) * pres_mod[34];
  //sp 88
  sp_rates[87] = -(fwd_rates[580] - rev_rates[534]) * pres_mod[34];

  //rxn 581
  //sp 61
  sp_rates[60] += (fwd_rates[581] - rev_rates[535]) * pres_mod[35];
  //sp 15
  sp_rates[14] += (fwd_rates[581] - rev_rates[535]) * pres_mod[35];
  //sp 88
  sp_rates[87] -= (fwd_rates[581] - rev_rates[535]) * pres_mod[35];

  //rxn 582
  //sp 89
  sp_rates[88] = (fwd_rates[582] - rev_rates[536]);
  //sp 4
  sp_rates[3] += (fwd_rates[582] - rev_rates[536]);
  //sp 88
  sp_rates[87] -= (fwd_rates[582] - rev_rates[536]);

  //rxn 583
  //sp 90
  sp_rates[89] = (fwd_rates[583] - rev_rates[537]);
  //sp 4
  sp_rates[3] += (fwd_rates[583] - rev_rates[537]);
  //sp 88
  sp_rates[87] -= (fwd_rates[583] - rev_rates[537]);

  //rxn 584
  //sp 89
  sp_rates[88] += (fwd_rates[584] - rev_rates[538]);
  //sp 3
  sp_rates[2] += (fwd_rates[584] - rev_rates[538]);
  //sp 4
  sp_rates[3] -= (fwd_rates[584] - rev_rates[538]);
  //sp 88
  sp_rates[87] -= (fwd_rates[584] - rev_rates[538]);

  //rxn 585
  //sp 90
  sp_rates[89] += (fwd_rates[585] - rev_rates[539]);
  //sp 3
  sp_rates[2] += (fwd_rates[585] - rev_rates[539]);
  //sp 4
  sp_rates[3] -= (fwd_rates[585] - rev_rates[539]);
  //sp 88
  sp_rates[87] -= (fwd_rates[585] - rev_rates[539]);

  //rxn 586
  //sp 89
  sp_rates[88] += (fwd_rates[586] - rev_rates[540]);
  //sp 7
  sp_rates[6] += (fwd_rates[586] - rev_rates[540]);
  //sp 88
  sp_rates[87] -= (fwd_rates[586] - rev_rates[540]);
  //sp 8
  sp_rates[7] -= (fwd_rates[586] - rev_rates[540]);

  //rxn 587
  //sp 90
  sp_rates[89] += (fwd_rates[587] - rev_rates[541]);
  //sp 7
  sp_rates[6] += (fwd_rates[587] - rev_rates[541]);
  //sp 88
  sp_rates[87] -= (fwd_rates[587] - rev_rates[541]);
  //sp 8
  sp_rates[7] -= (fwd_rates[587] - rev_rates[541]);

  //rxn 588
  //sp 10
  sp_rates[9] += (fwd_rates[588] - rev_rates[542]);
  //sp 5
  sp_rates[4] -= (fwd_rates[588] - rev_rates[542]);
  //sp 62
  sp_rates[61] += (fwd_rates[588] - rev_rates[542]);
  //sp 61
  sp_rates[60] -= (fwd_rates[588] - rev_rates[542]);

  //rxn 589
  //sp 89
  sp_rates[88] -= (fwd_rates[589] - rev_rates[543]);
  //sp 10
  sp_rates[9] += (fwd_rates[589] - rev_rates[543]);
  //sp 5
  sp_rates[4] -= (fwd_rates[589] - rev_rates[543]);
  //sp 97
  sp_rates[96] = (fwd_rates[589] - rev_rates[543]);

  //rxn 590
  //sp 97
  sp_rates[96] += (fwd_rates[590] - rev_rates[544]);
  //sp 90
  sp_rates[89] -= (fwd_rates[590] - rev_rates[544]);
  //sp 5
  sp_rates[4] -= (fwd_rates[590] - rev_rates[544]);
  //sp 10
  sp_rates[9] += (fwd_rates[590] - rev_rates[544]);

  //rxn 591
  //sp 90
  sp_rates[89] -= (fwd_rates[591] - rev_rates[545]);
  //sp 5
  sp_rates[4] -= (fwd_rates[591] - rev_rates[545]);
  //sp 10
  sp_rates[9] += (fwd_rates[591] - rev_rates[545]);
  //sp 98
  sp_rates[97] = (fwd_rates[591] - rev_rates[545]);

  //rxn 592
  //sp 91
  sp_rates[90] = -(fwd_rates[592] - rev_rates[546]) * pres_mod[36];
  //sp 60
  sp_rates[59] += (fwd_rates[592] - rev_rates[546]) * pres_mod[36];
  //sp 15
  sp_rates[14] += (fwd_rates[592] - rev_rates[546]) * pres_mod[36];

  //rxn 593
  //sp 91
  sp_rates[90] -= (fwd_rates[593] - rev_rates[547]);
  //sp 4
  sp_rates[3] += (fwd_rates[593] - rev_rates[547]);
  //sp 93
  sp_rates[92] = (fwd_rates[593] - rev_rates[547]);

  //rxn 594
  //sp 92
  sp_rates[91] = (fwd_rates[594] - rev_rates[548]);
  //sp 91
  sp_rates[90] -= (fwd_rates[594] - rev_rates[548]);
  //sp 4
  sp_rates[3] += (fwd_rates[594] - rev_rates[548]);

  //rxn 595
  //sp 91
  sp_rates[90] -= (fwd_rates[595] - rev_rates[549]);
  //sp 92
  sp_rates[91] += (fwd_rates[595] - rev_rates[549]);
  //sp 14
  sp_rates[13] += (fwd_rates[595] - rev_rates[549]);
  //sp 15
  sp_rates[14] -= (fwd_rates[595] - rev_rates[549]);

  //rxn 596
  //sp 3
  sp_rates[2] += (fwd_rates[596] - rev_rates[550]);
  //sp 92
  sp_rates[91] += (fwd_rates[596] - rev_rates[550]);
  //sp 91
  sp_rates[90] -= (fwd_rates[596] - rev_rates[550]);
  //sp 4
  sp_rates[3] -= (fwd_rates[596] - rev_rates[550]);

  //rxn 597
  //sp 91
  sp_rates[90] -= (fwd_rates[597] - rev_rates[551]);
  //sp 92
  sp_rates[91] += (fwd_rates[597] - rev_rates[551]);
  //sp 7
  sp_rates[6] += (fwd_rates[597] - rev_rates[551]);
  //sp 8
  sp_rates[7] -= (fwd_rates[597] - rev_rates[551]);

  //rxn 598
  //sp 91
  sp_rates[90] -= (fwd_rates[598] - rev_rates[552]);
  //sp 92
  sp_rates[91] += (fwd_rates[598] - rev_rates[552]);
  //sp 30
  sp_rates[29] += (fwd_rates[598] - rev_rates[552]);
  //sp 31
  sp_rates[30] -= (fwd_rates[598] - rev_rates[552]);

  //rxn 599
  //sp 9
  sp_rates[8] += (fwd_rates[599] - rev_rates[553]);
  //sp 10
  sp_rates[9] -= (fwd_rates[599] - rev_rates[553]);
  //sp 91
  sp_rates[90] -= (fwd_rates[599] - rev_rates[553]);
  //sp 92
  sp_rates[91] += (fwd_rates[599] - rev_rates[553]);

  //rxn 600
  //sp 91
  sp_rates[90] -= (fwd_rates[600] - rev_rates[554]);
  //sp 92
  sp_rates[91] += (fwd_rates[600] - rev_rates[554]);
  //sp 6
  sp_rates[5] -= (fwd_rates[600] - rev_rates[554]);
  //sp 8
  sp_rates[7] += (fwd_rates[600] - rev_rates[554]);

  //rxn 601
  //sp 91
  sp_rates[90] -= (fwd_rates[601] - rev_rates[555]);
  //sp 92
  sp_rates[91] += (fwd_rates[601] - rev_rates[555]);
  //sp 22
  sp_rates[21] += (fwd_rates[601] - rev_rates[555]);
  //sp 23
  sp_rates[22] -= (fwd_rates[601] - rev_rates[555]);

  //rxn 602
  //sp 10
  sp_rates[9] += (fwd_rates[602] - rev_rates[556]);
  //sp 91
  sp_rates[90] -= (fwd_rates[602] - rev_rates[556]);
  //sp 92
  sp_rates[91] += (fwd_rates[602] - rev_rates[556]);
  //sp 5
  sp_rates[4] -= (fwd_rates[602] - rev_rates[556]);

  //rxn 603
  //sp 92
  sp_rates[91] += (fwd_rates[603] - rev_rates[557]);
  //sp 91
  sp_rates[90] -= (fwd_rates[603] - rev_rates[557]);
  //sp 20
  sp_rates[19] += (fwd_rates[603] - rev_rates[557]);
  //sp 21
  sp_rates[20] -= (fwd_rates[603] - rev_rates[557]);

  //rxn 604
  //sp 33
  sp_rates[32] -= (fwd_rates[604] - rev_rates[558]);
  //sp 91
  sp_rates[90] -= (fwd_rates[604] - rev_rates[558]);
  //sp 92
  sp_rates[91] += (fwd_rates[604] - rev_rates[558]);
  //sp 32
  sp_rates[31] += (fwd_rates[604] - rev_rates[558]);

  //rxn 605
  //sp 3
  sp_rates[2] += (fwd_rates[605] - rev_rates[559]);
  //sp 91
  sp_rates[90] -= (fwd_rates[605] - rev_rates[559]);
  //sp 4
  sp_rates[3] -= (fwd_rates[605] - rev_rates[559]);
  //sp 93
  sp_rates[92] += (fwd_rates[605] - rev_rates[559]);

  //rxn 606
  //sp 91
  sp_rates[90] -= (fwd_rates[606] - rev_rates[560]);
  //sp 93
  sp_rates[92] += (fwd_rates[606] - rev_rates[560]);
  //sp 14
  sp_rates[13] += (fwd_rates[606] - rev_rates[560]);
  //sp 15
  sp_rates[14] -= (fwd_rates[606] - rev_rates[560]);

  //rxn 607
  //sp 91
  sp_rates[90] -= (fwd_rates[607] - rev_rates[561]);
  //sp 93
  sp_rates[92] += (fwd_rates[607] - rev_rates[561]);
  //sp 7
  sp_rates[6] += (fwd_rates[607] - rev_rates[561]);
  //sp 8
  sp_rates[7] -= (fwd_rates[607] - rev_rates[561]);

  //rxn 608
  //sp 91
  sp_rates[90] -= (fwd_rates[608] - rev_rates[562]);
  //sp 93
  sp_rates[92] += (fwd_rates[608] - rev_rates[562]);
  //sp 30
  sp_rates[29] += (fwd_rates[608] - rev_rates[562]);
  //sp 31
  sp_rates[30] -= (fwd_rates[608] - rev_rates[562]);

  //rxn 609
  //sp 9
  sp_rates[8] += (fwd_rates[609] - rev_rates[563]);
  //sp 10
  sp_rates[9] -= (fwd_rates[609] - rev_rates[563]);
  //sp 91
  sp_rates[90] -= (fwd_rates[609] - rev_rates[563]);
  //sp 93
  sp_rates[92] += (fwd_rates[609] - rev_rates[563]);

  //rxn 610
  //sp 91
  sp_rates[90] -= (fwd_rates[610] - rev_rates[564]);
  //sp 93
  sp_rates[92] += (fwd_rates[610] - rev_rates[564]);
  //sp 6
  sp_rates[5] -= (fwd_rates[610] - rev_rates[564]);
  //sp 8
  sp_rates[7] += (fwd_rates[610] - rev_rates[564]);

  //rxn 611
  //sp 91
  sp_rates[90] -= (fwd_rates[611] - rev_rates[565]);
  //sp 93
  sp_rates[92] += (fwd_rates[611] - rev_rates[565]);
  //sp 22
  sp_rates[21] += (fwd_rates[611] - rev_rates[565]);
  //sp 23
  sp_rates[22] -= (fwd_rates[611] - rev_rates[565]);

  //rxn 612
  //sp 10
  sp_rates[9] += (fwd_rates[612] - rev_rates[566]);
  //sp 91
  sp_rates[90] -= (fwd_rates[612] - rev_rates[566]);
  //sp 5
  sp_rates[4] -= (fwd_rates[612] - rev_rates[566]);
  //sp 93
  sp_rates[92] += (fwd_rates[612] - rev_rates[566]);

  //rxn 613
  //sp 91
  sp_rates[90] -= (fwd_rates[613] - rev_rates[567]);
  //sp 20
  sp_rates[19] += (fwd_rates[613] - rev_rates[567]);
  //sp 21
  sp_rates[20] -= (fwd_rates[613] - rev_rates[567]);
  //sp 93
  sp_rates[92] += (fwd_rates[613] - rev_rates[567]);

  //rxn 614
  //sp 33
  sp_rates[32] -= (fwd_rates[614] - rev_rates[568]);
  //sp 91
  sp_rates[90] -= (fwd_rates[614] - rev_rates[568]);
  //sp 93
  sp_rates[92] += (fwd_rates[614] - rev_rates[568]);
  //sp 32
  sp_rates[31] += (fwd_rates[614] - rev_rates[568]);

  //rxn 615
  //sp 92
  sp_rates[91] -= (fwd_rates[615] - rev_rates[569]);
  //sp 93
  sp_rates[92] += (fwd_rates[615] - rev_rates[569]);

  //rxn 616
  //sp 92
  sp_rates[91] -= (fwd_rates[616] - rev_rates[570]);
  //sp 93
  sp_rates[92] += (fwd_rates[616] - rev_rates[570]);

  //rxn 617
  //sp 92
  sp_rates[91] -= (fwd_rates[617] - rev_rates[571]);
  //sp 62
  sp_rates[61] += (fwd_rates[617] - rev_rates[571]);
  //sp 15
  sp_rates[14] += (fwd_rates[617] - rev_rates[571]);

  //rxn 618
  //sp 65
  sp_rates[64] += (fwd_rates[618] - rev_rates[572]);
  //sp 94
  sp_rates[93] = -(fwd_rates[618] - rev_rates[572]);
  //sp 15
  sp_rates[14] += (fwd_rates[618] - rev_rates[572]);

  //rxn 619
  //sp 4
  sp_rates[3] += (fwd_rates[619] - rev_rates[573]);
  //sp 94
  sp_rates[93] -= (fwd_rates[619] - rev_rates[573]);
  //sp 95
  sp_rates[94] = (fwd_rates[619] - rev_rates[573]);

  //rxn 620
  //sp 95
  sp_rates[94] += (fwd_rates[620] - rev_rates[574]);
  //sp 94
  sp_rates[93] -= (fwd_rates[620] - rev_rates[574]);
  //sp 7
  sp_rates[6] += (fwd_rates[620] - rev_rates[574]);
  //sp 8
  sp_rates[7] -= (fwd_rates[620] - rev_rates[574]);

  //rxn 621
  //sp 10
  sp_rates[9] += (fwd_rates[621] - rev_rates[575]);
  //sp 5
  sp_rates[4] -= (fwd_rates[621] - rev_rates[575]);
  //sp 94
  sp_rates[93] -= (fwd_rates[621] - rev_rates[575]);
  //sp 95
  sp_rates[94] += (fwd_rates[621] - rev_rates[575]);

  //rxn 622
  //sp 94
  sp_rates[93] -= (fwd_rates[622] - rev_rates[576]);
  //sp 6
  sp_rates[5] -= (fwd_rates[622] - rev_rates[576]);
  //sp 95
  sp_rates[94] += (fwd_rates[622] - rev_rates[576]);
  //sp 8
  sp_rates[7] += (fwd_rates[622] - rev_rates[576]);

  //rxn 623
  //sp 9
  sp_rates[8] += (fwd_rates[623] - rev_rates[577]);
  //sp 10
  sp_rates[9] -= (fwd_rates[623] - rev_rates[577]);
  //sp 94
  sp_rates[93] -= (fwd_rates[623] - rev_rates[577]);
  //sp 95
  sp_rates[94] += (fwd_rates[623] - rev_rates[577]);

  //rxn 624
  //sp 95
  sp_rates[94] += (fwd_rates[624] - rev_rates[578]);
  //sp 22
  sp_rates[21] += (fwd_rates[624] - rev_rates[578]);
  //sp 94
  sp_rates[93] -= (fwd_rates[624] - rev_rates[578]);
  //sp 23
  sp_rates[22] -= (fwd_rates[624] - rev_rates[578]);

  //rxn 625
  //sp 95
  sp_rates[94] += (fwd_rates[625] - rev_rates[579]);
  //sp 62
  sp_rates[61] += (fwd_rates[625] - rev_rates[579]);
  //sp 94
  sp_rates[93] -= (fwd_rates[625] - rev_rates[579]);
  //sp 63
  sp_rates[62] -= (fwd_rates[625] - rev_rates[579]);

  //rxn 626
  //sp 62
  sp_rates[61] += (fwd_rates[626] - rev_rates[580]);
  //sp 94
  sp_rates[93] -= (fwd_rates[626] - rev_rates[580]);
  //sp 95
  sp_rates[94] += (fwd_rates[626] - rev_rates[580]);
  //sp 64
  sp_rates[63] -= (fwd_rates[626] - rev_rates[580]);

  //rxn 627
  //sp 65
  sp_rates[64] -= (fwd_rates[627] - rev_rates[581]);
  //sp 62
  sp_rates[61] += (fwd_rates[627] - rev_rates[581]);
  //sp 94
  sp_rates[93] -= (fwd_rates[627] - rev_rates[581]);
  //sp 95
  sp_rates[94] += (fwd_rates[627] - rev_rates[581]);

  //rxn 628
  //sp 20
  sp_rates[19] += (fwd_rates[628] - rev_rates[582]);
  //sp 21
  sp_rates[20] -= (fwd_rates[628] - rev_rates[582]);
  //sp 94
  sp_rates[93] -= (fwd_rates[628] - rev_rates[582]);
  //sp 95
  sp_rates[94] += (fwd_rates[628] - rev_rates[582]);

  //rxn 629
  //sp 79
  sp_rates[78] -= (fwd_rates[629] - rev_rates[583]);
  //sp 15
  sp_rates[14] -= (fwd_rates[629] - rev_rates[583]);
  //sp 95
  sp_rates[94] += (fwd_rates[629] - rev_rates[583]);

  //rxn 630
  //sp 10
  sp_rates[9] -= (fwd_rates[630] - rev_rates[584]);
  //sp 96
  sp_rates[95] = (fwd_rates[630] - rev_rates[584]);
  //sp 95
  sp_rates[94] -= (fwd_rates[630] - rev_rates[584]);
  //sp 8
  sp_rates[7] += (fwd_rates[630] - rev_rates[584]);

  //rxn 631
  //sp 25
  sp_rates[24] += (fwd_rates[631] - rev_rates[585]);
  //sp 65
  sp_rates[64] += (fwd_rates[631] - rev_rates[585]);
  //sp 96
  sp_rates[95] -= (fwd_rates[631] - rev_rates[585]);

  //rxn 632
  //sp 26
  sp_rates[25] += (fwd_rates[632] - rev_rates[586]);
  //sp 62
  sp_rates[61] += (fwd_rates[632] - rev_rates[586]);
  //sp 96
  sp_rates[95] -= (fwd_rates[632] - rev_rates[586]);

  //rxn 633
  //sp 25
  sp_rates[24] -= (fwd_rates[633] - rev_rates[587]);
  //sp 26
  sp_rates[25] += (fwd_rates[633] - rev_rates[587]);
  //sp 65
  sp_rates[64] -= (fwd_rates[633] - rev_rates[587]);
  //sp 62
  sp_rates[61] += (fwd_rates[633] - rev_rates[587]);

  //rxn 634
  //sp 21
  sp_rates[20] -= (fwd_rates[634] - rev_rates[588]);
  //sp 96
  sp_rates[95] += (fwd_rates[634] - rev_rates[588]);
  //sp 95
  sp_rates[94] -= (fwd_rates[634] - rev_rates[588]);
  //sp 23
  sp_rates[22] += (fwd_rates[634] - rev_rates[588]);

  //rxn 635
  //sp 92
  sp_rates[91] += (fwd_rates[635] - rev_rates[589]);
  //sp 4
  sp_rates[3] -= (fwd_rates[635] - rev_rates[589]);
  //sp 94
  sp_rates[93] -= (fwd_rates[635] - rev_rates[589]);

  //rxn 636
  //sp 92
  sp_rates[91] += (fwd_rates[636] - rev_rates[590]);
  //sp 4
  sp_rates[3] -= (fwd_rates[636] - rev_rates[590]);
  //sp 94
  sp_rates[93] -= (fwd_rates[636] - rev_rates[590]);

  //rxn 637
  //sp 4
  sp_rates[3] -= (fwd_rates[637] - rev_rates[591]);
  //sp 93
  sp_rates[92] += (fwd_rates[637] - rev_rates[591]);
  //sp 94
  sp_rates[93] -= (fwd_rates[637] - rev_rates[591]);

  //rxn 638
  //sp 4
  sp_rates[3] -= (fwd_rates[638] - rev_rates[592]);
  //sp 93
  sp_rates[92] += (fwd_rates[638] - rev_rates[592]);
  //sp 94
  sp_rates[93] -= (fwd_rates[638] - rev_rates[592]);

  //rxn 639
  //sp 62
  sp_rates[61] += (fwd_rates[639] - rev_rates[593]);
  //sp 4
  sp_rates[3] -= (fwd_rates[639] - rev_rates[593]);
  //sp 94
  sp_rates[93] -= (fwd_rates[639] - rev_rates[593]);
  //sp 15
  sp_rates[14] += (fwd_rates[639] - rev_rates[593]);

  //rxn 640
  //sp 62
  sp_rates[61] += (fwd_rates[640] - rev_rates[594]);
  //sp 4
  sp_rates[3] -= (fwd_rates[640] - rev_rates[594]);
  //sp 94
  sp_rates[93] -= (fwd_rates[640] - rev_rates[594]);
  //sp 15
  sp_rates[14] += (fwd_rates[640] - rev_rates[594]);

  //rxn 641
  //sp 86
  sp_rates[85] += (fwd_rates[641] - rev_rates[595]);
  //sp 94
  sp_rates[93] -= (fwd_rates[641] - rev_rates[595]);
  //sp 15
  sp_rates[14] += (fwd_rates[641] - rev_rates[595]);
  //sp 8
  sp_rates[7] -= (fwd_rates[641] - rev_rates[595]);

  //rxn 642
  //sp 63
  sp_rates[62] -= (fwd_rates[642] - rev_rates[596]) * pres_mod[37];
  //sp 97
  sp_rates[96] += (fwd_rates[642] - rev_rates[596]) * pres_mod[37];
  //sp 15
  sp_rates[14] -= (fwd_rates[642] - rev_rates[596]) * pres_mod[37];

  //rxn 643
  //sp 97
  sp_rates[96] += (fwd_rates[643] - rev_rates[597]) * pres_mod[38];
  //sp 35
  sp_rates[34] -= (fwd_rates[643] - rev_rates[597]) * pres_mod[38];
  //sp 31
  sp_rates[30] -= (fwd_rates[643] - rev_rates[597]) * pres_mod[38];

  //rxn 644
  //sp 100
  sp_rates[99] = -(fwd_rates[644] - rev_rates[598]) * pres_mod[39];
  //sp 4
  sp_rates[3] -= (fwd_rates[644] - rev_rates[598]) * pres_mod[39];
  //sp 97
  sp_rates[96] += (fwd_rates[644] - rev_rates[598]) * pres_mod[39];

  //rxn 645
  //sp 97
  sp_rates[96] += (fwd_rates[645] - rev_rates[599]);
  //sp 99
  sp_rates[98] = -(fwd_rates[645] - rev_rates[599]);
  //sp 4
  sp_rates[3] -= (fwd_rates[645] - rev_rates[599]);

  //rxn 646
  //sp 98
  sp_rates[97] += (fwd_rates[646] - rev_rates[600]) * pres_mod[40];
  //sp 99
  sp_rates[98] -= (fwd_rates[646] - rev_rates[600]) * pres_mod[40];
  //sp 4
  sp_rates[3] -= (fwd_rates[646] - rev_rates[600]) * pres_mod[40];

  //rxn 647
  //sp 98
  sp_rates[97] += (fwd_rates[647] - rev_rates[601]) * pres_mod[41];
  //sp 15
  sp_rates[14] -= (fwd_rates[647] - rev_rates[601]) * pres_mod[41];
  //sp 64
  sp_rates[63] -= (fwd_rates[647] - rev_rates[601]) * pres_mod[41];

  //rxn 648
  //sp 63
  sp_rates[62] += (fwd_rates[648] - rev_rates[602]);
  //sp 98
  sp_rates[97] -= (fwd_rates[648] - rev_rates[602]);
  //sp 15
  sp_rates[14] += (fwd_rates[648] - rev_rates[602]);

  //rxn 649
  //sp 98
  sp_rates[97] -= (fwd_rates[649] - rev_rates[603]);
  //sp 99
  sp_rates[98] += (fwd_rates[649] - rev_rates[603]);
  //sp 4
  sp_rates[3] += (fwd_rates[649] - rev_rates[603]);

  //rxn 650
  //sp 97
  sp_rates[96] -= (fwd_rates[650] - rev_rates[604]);
  //sp 99
  sp_rates[98] += (fwd_rates[650] - rev_rates[604]);
  //sp 7
  sp_rates[6] += (fwd_rates[650] - rev_rates[604]);
  //sp 8
  sp_rates[7] -= (fwd_rates[650] - rev_rates[604]);

  //rxn 651
  //sp 97
  sp_rates[96] -= (fwd_rates[651] - rev_rates[605]);
  //sp 100
  sp_rates[99] += (fwd_rates[651] - rev_rates[605]);
  //sp 7
  sp_rates[6] += (fwd_rates[651] - rev_rates[605]);
  //sp 8
  sp_rates[7] -= (fwd_rates[651] - rev_rates[605]);

  //rxn 652
  //sp 98
  sp_rates[97] -= (fwd_rates[652] - rev_rates[606]);
  //sp 99
  sp_rates[98] += (fwd_rates[652] - rev_rates[606]);
  //sp 7
  sp_rates[6] += (fwd_rates[652] - rev_rates[606]);
  //sp 8
  sp_rates[7] -= (fwd_rates[652] - rev_rates[606]);

  //rxn 653
  //sp 97
  sp_rates[96] -= (fwd_rates[653] - rev_rates[607]);
  //sp 10
  sp_rates[9] -= (fwd_rates[653] - rev_rates[607]);
  //sp 99
  sp_rates[98] += (fwd_rates[653] - rev_rates[607]);
  //sp 9
  sp_rates[8] += (fwd_rates[653] - rev_rates[607]);

  //rxn 654
  //sp 9
  sp_rates[8] += (fwd_rates[654] - rev_rates[608]);
  //sp 10
  sp_rates[9] -= (fwd_rates[654] - rev_rates[608]);
  //sp 99
  sp_rates[98] += (fwd_rates[654] - rev_rates[608]);
  //sp 98
  sp_rates[97] -= (fwd_rates[654] - rev_rates[608]);

  //rxn 655
  //sp 97
  sp_rates[96] -= (fwd_rates[655] - rev_rates[609]);
  //sp 10
  sp_rates[9] -= (fwd_rates[655] - rev_rates[609]);
  //sp 100
  sp_rates[99] += (fwd_rates[655] - rev_rates[609]);
  //sp 9
  sp_rates[8] += (fwd_rates[655] - rev_rates[609]);

  //rxn 656
  //sp 98
  sp_rates[97] -= (fwd_rates[656] - rev_rates[610]);
  //sp 3
  sp_rates[2] += (fwd_rates[656] - rev_rates[610]);
  //sp 4
  sp_rates[3] -= (fwd_rates[656] - rev_rates[610]);
  //sp 99
  sp_rates[98] += (fwd_rates[656] - rev_rates[610]);

  //rxn 657
  //sp 97
  sp_rates[96] -= (fwd_rates[657] - rev_rates[611]);
  //sp 99
  sp_rates[98] += (fwd_rates[657] - rev_rates[611]);
  //sp 6
  sp_rates[5] -= (fwd_rates[657] - rev_rates[611]);
  //sp 8
  sp_rates[7] += (fwd_rates[657] - rev_rates[611]);

  //rxn 658
  //sp 98
  sp_rates[97] -= (fwd_rates[658] - rev_rates[612]);
  //sp 99
  sp_rates[98] += (fwd_rates[658] - rev_rates[612]);
  //sp 6
  sp_rates[5] -= (fwd_rates[658] - rev_rates[612]);
  //sp 8
  sp_rates[7] += (fwd_rates[658] - rev_rates[612]);

  //rxn 659
  //sp 97
  sp_rates[96] -= (fwd_rates[659] - rev_rates[613]);
  //sp 100
  sp_rates[99] += (fwd_rates[659] - rev_rates[613]);
  //sp 6
  sp_rates[5] -= (fwd_rates[659] - rev_rates[613]);
  //sp 8
  sp_rates[7] += (fwd_rates[659] - rev_rates[613]);

  //rxn 660
  //sp 97
  sp_rates[96] -= (fwd_rates[660] - rev_rates[614]);
  //sp 10
  sp_rates[9] += (fwd_rates[660] - rev_rates[614]);
  //sp 99
  sp_rates[98] += (fwd_rates[660] - rev_rates[614]);
  //sp 5
  sp_rates[4] -= (fwd_rates[660] - rev_rates[614]);

  //rxn 661
  //sp 98
  sp_rates[97] -= (fwd_rates[661] - rev_rates[615]);
  //sp 99
  sp_rates[98] += (fwd_rates[661] - rev_rates[615]);
  //sp 5
  sp_rates[4] -= (fwd_rates[661] - rev_rates[615]);
  //sp 10
  sp_rates[9] += (fwd_rates[661] - rev_rates[615]);

  //rxn 662
  //sp 97
  sp_rates[96] -= (fwd_rates[662] - rev_rates[616]);
  //sp 10
  sp_rates[9] += (fwd_rates[662] - rev_rates[616]);
  //sp 100
  sp_rates[99] += (fwd_rates[662] - rev_rates[616]);
  //sp 5
  sp_rates[4] -= (fwd_rates[662] - rev_rates[616]);

  //rxn 663
  //sp 98
  sp_rates[97] -= (fwd_rates[663] - rev_rates[617]);
  //sp 99
  sp_rates[98] += (fwd_rates[663] - rev_rates[617]);
  //sp 14
  sp_rates[13] += (fwd_rates[663] - rev_rates[617]);
  //sp 15
  sp_rates[14] -= (fwd_rates[663] - rev_rates[617]);

  //rxn 664
  //sp 97
  sp_rates[96] -= (fwd_rates[664] - rev_rates[618]);
  //sp 99
  sp_rates[98] += (fwd_rates[664] - rev_rates[618]);
  //sp 20
  sp_rates[19] += (fwd_rates[664] - rev_rates[618]);
  //sp 21
  sp_rates[20] -= (fwd_rates[664] - rev_rates[618]);

  //rxn 665
  //sp 97
  sp_rates[96] -= (fwd_rates[665] - rev_rates[619]);
  //sp 100
  sp_rates[99] += (fwd_rates[665] - rev_rates[619]);
  //sp 20
  sp_rates[19] += (fwd_rates[665] - rev_rates[619]);
  //sp 21
  sp_rates[20] -= (fwd_rates[665] - rev_rates[619]);

  //rxn 666
  //sp 97
  sp_rates[96] -= (fwd_rates[666] - rev_rates[620]);
  //sp 99
  sp_rates[98] += (fwd_rates[666] - rev_rates[620]);
  //sp 22
  sp_rates[21] += (fwd_rates[666] - rev_rates[620]);
  //sp 23
  sp_rates[22] -= (fwd_rates[666] - rev_rates[620]);

  //rxn 667
  //sp 97
  sp_rates[96] -= (fwd_rates[667] - rev_rates[621]);
  //sp 100
  sp_rates[99] += (fwd_rates[667] - rev_rates[621]);
  //sp 22
  sp_rates[21] += (fwd_rates[667] - rev_rates[621]);
  //sp 23
  sp_rates[22] -= (fwd_rates[667] - rev_rates[621]);

  //rxn 668
  //sp 98
  sp_rates[97] -= (fwd_rates[668] - rev_rates[622]);
  //sp 99
  sp_rates[98] += (fwd_rates[668] - rev_rates[622]);
  //sp 22
  sp_rates[21] += (fwd_rates[668] - rev_rates[622]);
  //sp 23
  sp_rates[22] -= (fwd_rates[668] - rev_rates[622]);

  //rxn 669
  //sp 33
  sp_rates[32] -= (fwd_rates[669] - rev_rates[623]);
  //sp 99
  sp_rates[98] += (fwd_rates[669] - rev_rates[623]);
  //sp 97
  sp_rates[96] -= (fwd_rates[669] - rev_rates[623]);
  //sp 32
  sp_rates[31] += (fwd_rates[669] - rev_rates[623]);

  //rxn 670
  //sp 97
  sp_rates[96] += (fwd_rates[670] - rev_rates[624]);
  //sp 34
  sp_rates[33] += (fwd_rates[670] - rev_rates[624]);
  //sp 99
  sp_rates[98] -= (fwd_rates[670] - rev_rates[624]);
  //sp 31
  sp_rates[30] -= (fwd_rates[670] - rev_rates[624]);

  //rxn 671
  //sp 25
  sp_rates[24] += (fwd_rates[671] - rev_rates[625]);
  //sp 99
  sp_rates[98] -= (fwd_rates[671] - rev_rates[625]);
  //sp 97
  sp_rates[96] += (fwd_rates[671] - rev_rates[625]);
  //sp 23
  sp_rates[22] -= (fwd_rates[671] - rev_rates[625]);

  //rxn 672
  //sp 99
  sp_rates[98] -= (fwd_rates[672] - rev_rates[626]);
  //sp 4
  sp_rates[3] += (fwd_rates[672] - rev_rates[626]);
  //sp 102
  sp_rates[101] = (fwd_rates[672] - rev_rates[626]);

  //rxn 673
  //sp 34
  sp_rates[33] += (fwd_rates[673] - rev_rates[627]);
  //sp 35
  sp_rates[34] += (fwd_rates[673] - rev_rates[627]);
  //sp 100
  sp_rates[99] -= (fwd_rates[673] - rev_rates[627]);

  //rxn 674
  //sp 4
  sp_rates[3] += (fwd_rates[674] - rev_rates[628]);
  //sp 100
  sp_rates[99] -= (fwd_rates[674] - rev_rates[628]);
  //sp 102
  sp_rates[101] += (fwd_rates[674] - rev_rates[628]);

  //rxn 675
  //sp 99
  sp_rates[98] -= (fwd_rates[675] - rev_rates[629]);
  //sp 100
  sp_rates[99] += (fwd_rates[675] - rev_rates[629]);

  //rxn 676
  //sp 10
  sp_rates[9] -= (fwd_rates[676] - rev_rates[630]);
  //sp 99
  sp_rates[98] -= (fwd_rates[676] - rev_rates[630]);
  //sp 101
  sp_rates[100] = (fwd_rates[676] - rev_rates[630]);
  //sp 8
  sp_rates[7] += (fwd_rates[676] - rev_rates[630]);

  //rxn 677
  //sp 25
  sp_rates[24] -= (fwd_rates[677] - rev_rates[631]);
  //sp 26
  sp_rates[25] += (fwd_rates[677] - rev_rates[631]);
  //sp 62
  sp_rates[61] += (fwd_rates[677] - rev_rates[631]);
  //sp 64
  sp_rates[63] -= (fwd_rates[677] - rev_rates[631]);

  //rxn 678
  //sp 101
  sp_rates[100] += (fwd_rates[678] - rev_rates[632]);
  //sp 99
  sp_rates[98] -= (fwd_rates[678] - rev_rates[632]);
  //sp 21
  sp_rates[20] -= (fwd_rates[678] - rev_rates[632]);
  //sp 23
  sp_rates[22] += (fwd_rates[678] - rev_rates[632]);

  //rxn 679
  //sp 33
  sp_rates[32] -= (fwd_rates[679] - rev_rates[633]);
  //sp 99
  sp_rates[98] -= (fwd_rates[679] - rev_rates[633]);
  //sp 101
  sp_rates[100] += (fwd_rates[679] - rev_rates[633]);
  //sp 39
  sp_rates[38] += (fwd_rates[679] - rev_rates[633]);

  //rxn 680
  //sp 101
  sp_rates[100] -= (fwd_rates[680] - rev_rates[634]);
  //sp 15
  sp_rates[14] += (fwd_rates[680] - rev_rates[634]);
  //sp 85
  sp_rates[84] += (fwd_rates[680] - rev_rates[634]);

  //rxn 681
  //sp 99
  sp_rates[98] -= (fwd_rates[681] - rev_rates[635]);
  //sp 85
  sp_rates[84] += (fwd_rates[681] - rev_rates[635]);
  //sp 6
  sp_rates[5] -= (fwd_rates[681] - rev_rates[635]);
  //sp 15
  sp_rates[14] += (fwd_rates[681] - rev_rates[635]);

  //rxn 682
  //sp 97
  sp_rates[96] -= (fwd_rates[682] - rev_rates[636]);
  //sp 34
  sp_rates[33] += (fwd_rates[682] - rev_rates[636]);
  //sp 4
  sp_rates[3] -= (fwd_rates[682] - rev_rates[636]);
  //sp 31
  sp_rates[30] += (fwd_rates[682] - rev_rates[636]);

  //rxn 683
  //sp 97
  sp_rates[96] -= (fwd_rates[683] - rev_rates[637]);
  //sp 34
  sp_rates[33] += (fwd_rates[683] - rev_rates[637]);
  //sp 4
  sp_rates[3] -= (fwd_rates[683] - rev_rates[637]);
  //sp 31
  sp_rates[30] += (fwd_rates[683] - rev_rates[637]);

  //rxn 684
  //sp 97
  sp_rates[96] -= (fwd_rates[684] - rev_rates[638]);
  //sp 4
  sp_rates[3] -= (fwd_rates[684] - rev_rates[638]);
  //sp 62
  sp_rates[61] += (fwd_rates[684] - rev_rates[638]);
  //sp 15
  sp_rates[14] += (fwd_rates[684] - rev_rates[638]);

  //rxn 685
  //sp 97
  sp_rates[96] -= (fwd_rates[685] - rev_rates[639]);
  //sp 4
  sp_rates[3] -= (fwd_rates[685] - rev_rates[639]);
  //sp 62
  sp_rates[61] += (fwd_rates[685] - rev_rates[639]);
  //sp 15
  sp_rates[14] += (fwd_rates[685] - rev_rates[639]);

  //rxn 686
  //sp 97
  sp_rates[96] -= (fwd_rates[686] - rev_rates[640]);
  //sp 89
  sp_rates[88] += (fwd_rates[686] - rev_rates[640]);
  //sp 4
  sp_rates[3] -= (fwd_rates[686] - rev_rates[640]);

  //rxn 687
  //sp 97
  sp_rates[96] -= (fwd_rates[687] - rev_rates[641]);
  //sp 89
  sp_rates[88] += (fwd_rates[687] - rev_rates[641]);
  //sp 4
  sp_rates[3] -= (fwd_rates[687] - rev_rates[641]);

  //rxn 688
  //sp 97
  sp_rates[96] -= (fwd_rates[688] - rev_rates[642]);
  //sp 90
  sp_rates[89] += (fwd_rates[688] - rev_rates[642]);
  //sp 4
  sp_rates[3] -= (fwd_rates[688] - rev_rates[642]);

  //rxn 689
  //sp 97
  sp_rates[96] -= (fwd_rates[689] - rev_rates[643]);
  //sp 90
  sp_rates[89] += (fwd_rates[689] - rev_rates[643]);
  //sp 4
  sp_rates[3] -= (fwd_rates[689] - rev_rates[643]);

  //rxn 690
  //sp 98
  sp_rates[97] -= (fwd_rates[690] - rev_rates[644]);
  //sp 4
  sp_rates[3] -= (fwd_rates[690] - rev_rates[644]);
  //sp 31
  sp_rates[30] += (fwd_rates[690] - rev_rates[644]);
  //sp 34
  sp_rates[33] += (fwd_rates[690] - rev_rates[644]);

  //rxn 691
  //sp 98
  sp_rates[97] -= (fwd_rates[691] - rev_rates[645]);
  //sp 4
  sp_rates[3] -= (fwd_rates[691] - rev_rates[645]);
  //sp 62
  sp_rates[61] += (fwd_rates[691] - rev_rates[645]);
  //sp 15
  sp_rates[14] += (fwd_rates[691] - rev_rates[645]);

  //rxn 692
  //sp 89
  sp_rates[88] += (fwd_rates[692] - rev_rates[646]);
  //sp 98
  sp_rates[97] -= (fwd_rates[692] - rev_rates[646]);
  //sp 4
  sp_rates[3] -= (fwd_rates[692] - rev_rates[646]);

  //rxn 693
  //sp 90
  sp_rates[89] += (fwd_rates[693] - rev_rates[647]);
  //sp 98
  sp_rates[97] -= (fwd_rates[693] - rev_rates[647]);
  //sp 4
  sp_rates[3] -= (fwd_rates[693] - rev_rates[647]);

  //rxn 694
  //sp 97
  sp_rates[96] -= (fwd_rates[694] - rev_rates[648]);
  //sp 98
  sp_rates[97] += (fwd_rates[694] - rev_rates[648]);

  //rxn 695
  //sp 97
  sp_rates[96] -= (fwd_rates[695] - rev_rates[649]);
  //sp 98
  sp_rates[97] += (fwd_rates[695] - rev_rates[649]);

  //rxn 696
  //sp 89
  sp_rates[88] += (fwd_rates[696] - rev_rates[650]);
  //sp 90
  sp_rates[89] -= (fwd_rates[696] - rev_rates[650]);

  //rxn 697
  //sp 89
  sp_rates[88] -= (fwd_rates[697] - rev_rates[651]);
  //sp 34
  sp_rates[33] += (fwd_rates[697] - rev_rates[651]);
  //sp 31
  sp_rates[30] += (fwd_rates[697] - rev_rates[651]);

  //rxn 698
  //sp 89
  sp_rates[88] -= (fwd_rates[698] - rev_rates[652]);
  //sp 62
  sp_rates[61] += (fwd_rates[698] - rev_rates[652]);
  //sp 15
  sp_rates[14] += (fwd_rates[698] - rev_rates[652]);

  //rxn 699
  //sp 90
  sp_rates[89] -= (fwd_rates[699] - rev_rates[653]);
  //sp 31
  sp_rates[30] += (fwd_rates[699] - rev_rates[653]);
  //sp 34
  sp_rates[33] += (fwd_rates[699] - rev_rates[653]);

  //rxn 700
  //sp 90
  sp_rates[89] -= (fwd_rates[700] - rev_rates[654]);
  //sp 62
  sp_rates[61] += (fwd_rates[700] - rev_rates[654]);
  //sp 15
  sp_rates[14] += (fwd_rates[700] - rev_rates[654]);

  //rxn 701
  //sp 97
  sp_rates[96] -= (fwd_rates[701] - rev_rates[655]);
  //sp 26
  sp_rates[25] += (fwd_rates[701] - rev_rates[655]);
  //sp 61
  sp_rates[60] += (fwd_rates[701] - rev_rates[655]);
  //sp 6
  sp_rates[5] -= (fwd_rates[701] - rev_rates[655]);

  //rxn 702
  //sp 97
  sp_rates[96] -= fwd_rates[702];
  //sp 4
  sp_rates[3] += fwd_rates[702];
  //sp 6
  sp_rates[5] -= fwd_rates[702];
  //sp 48
  sp_rates[47] += fwd_rates[702];
  //sp 31
  sp_rates[30] += fwd_rates[702];

  //rxn 703
  //sp 98
  sp_rates[97] -= fwd_rates[703];
  //sp 6
  sp_rates[5] -= fwd_rates[703];
  //sp 11
  sp_rates[10] += fwd_rates[703];
  //sp 15
  sp_rates[14] += fwd_rates[703];
  //sp 31
  sp_rates[30] += fwd_rates[703];

  //rxn 704
  //sp 98
  sp_rates[97] -= fwd_rates[704];
  //sp 4
  sp_rates[3] += fwd_rates[704];
  //sp 6
  sp_rates[5] -= fwd_rates[704];
  //sp 48
  sp_rates[47] += fwd_rates[704];
  //sp 31
  sp_rates[30] += fwd_rates[704];

  //rxn 705
  //sp 97
  sp_rates[96] -= (fwd_rates[705] - rev_rates[656]);
  //sp 98
  sp_rates[97] += (fwd_rates[705] - rev_rates[656]);

  //rxn 706
  //sp 97
  sp_rates[96] -= (fwd_rates[706] - rev_rates[657]);
  //sp 3
  sp_rates[2] += (fwd_rates[706] - rev_rates[657]);
  //sp 102
  sp_rates[101] += (fwd_rates[706] - rev_rates[657]);

  //rxn 707
  //sp 98
  sp_rates[97] -= (fwd_rates[707] - rev_rates[658]);
  //sp 3
  sp_rates[2] += (fwd_rates[707] - rev_rates[658]);
  //sp 102
  sp_rates[101] += (fwd_rates[707] - rev_rates[658]);

  //rxn 708
  //sp 3
  sp_rates[2] += (fwd_rates[708] - rev_rates[659]);
  //sp 102
  sp_rates[101] -= (fwd_rates[708] - rev_rates[659]);
  //sp 104
  sp_rates[103] = (fwd_rates[708] - rev_rates[659]);

  //rxn 709
  //sp 34
  sp_rates[33] += (fwd_rates[709] - rev_rates[660]);
  //sp 35
  sp_rates[34] += (fwd_rates[709] - rev_rates[660]);
  //sp 4
  sp_rates[3] -= (fwd_rates[709] - rev_rates[660]);
  //sp 102
  sp_rates[101] -= (fwd_rates[709] - rev_rates[660]);

  //rxn 710
  //sp 4
  sp_rates[3] -= (fwd_rates[710] - rev_rates[661]);
  //sp 78
  sp_rates[77] += (fwd_rates[710] - rev_rates[661]);
  //sp 102
  sp_rates[101] -= (fwd_rates[710] - rev_rates[661]);
  //sp 15
  sp_rates[14] += (fwd_rates[710] - rev_rates[661]);

  //rxn 711
  //sp 79
  sp_rates[78] += (fwd_rates[711] - rev_rates[662]);
  //sp 4
  sp_rates[3] -= (fwd_rates[711] - rev_rates[662]);
  //sp 102
  sp_rates[101] -= (fwd_rates[711] - rev_rates[662]);
  //sp 15
  sp_rates[14] += (fwd_rates[711] - rev_rates[662]);

  //rxn 712
  //sp 102
  sp_rates[101] -= (fwd_rates[712] - rev_rates[663]);
  //sp 36
  sp_rates[35] += (fwd_rates[712] - rev_rates[663]);
  //sp 43
  sp_rates[42] += (fwd_rates[712] - rev_rates[663]);
  //sp 6
  sp_rates[5] -= (fwd_rates[712] - rev_rates[663]);

  //rxn 713
  //sp 85
  sp_rates[84] += (fwd_rates[713] - rev_rates[664]);
  //sp 102
  sp_rates[101] -= (fwd_rates[713] - rev_rates[664]);
  //sp 15
  sp_rates[14] += (fwd_rates[713] - rev_rates[664]);
  //sp 8
  sp_rates[7] -= (fwd_rates[713] - rev_rates[664]);

  //rxn 714
  //sp 25
  sp_rates[24] += (fwd_rates[714] - rev_rates[665]);
  //sp 102
  sp_rates[101] -= (fwd_rates[714] - rev_rates[665]);
  //sp 63
  sp_rates[62] += (fwd_rates[714] - rev_rates[665]);
  //sp 8
  sp_rates[7] -= (fwd_rates[714] - rev_rates[665]);

  //rxn 715
  //sp 4
  sp_rates[3] += (fwd_rates[715] - rev_rates[666]);
  //sp 35
  sp_rates[34] -= (fwd_rates[715] - rev_rates[666]);
  //sp 36
  sp_rates[35] -= (fwd_rates[715] - rev_rates[666]);
  //sp 104
  sp_rates[103] += (fwd_rates[715] - rev_rates[666]);

  //rxn 716
  //sp 35
  sp_rates[34] -= 2.0 * (fwd_rates[716] - rev_rates[667]);
  //sp 102
  sp_rates[101] += (fwd_rates[716] - rev_rates[667]);

  //rxn 717
  //sp 4
  sp_rates[3] += (fwd_rates[717] - rev_rates[668]);
  //sp 104
  sp_rates[103] += (fwd_rates[717] - rev_rates[668]);
  //sp 80
  sp_rates[79] -= (fwd_rates[717] - rev_rates[668]);
  //sp 16
  sp_rates[15] -= (fwd_rates[717] - rev_rates[668]);

  //rxn 718
  //sp 26
  sp_rates[25] += (fwd_rates[718] - rev_rates[669]);
  //sp 6
  sp_rates[5] -= (fwd_rates[718] - rev_rates[669]);
  //sp 80
  sp_rates[79] += (fwd_rates[718] - rev_rates[669]);
  //sp 104
  sp_rates[103] -= (fwd_rates[718] - rev_rates[669]);

  //rxn 719
  //sp 25
  sp_rates[24] += (fwd_rates[719] - rev_rates[670]);
  //sp 104
  sp_rates[103] -= (fwd_rates[719] - rev_rates[670]);
  //sp 80
  sp_rates[79] += (fwd_rates[719] - rev_rates[670]);
  //sp 8
  sp_rates[7] -= (fwd_rates[719] - rev_rates[670]);

  //rxn 720
  //sp 49
  sp_rates[48] -= (fwd_rates[720] - rev_rates[671]);
  //sp 11
  sp_rates[10] += (fwd_rates[720] - rev_rates[671]);
  //sp 104
  sp_rates[103] += (fwd_rates[720] - rev_rates[671]);
  //sp 80
  sp_rates[79] -= (fwd_rates[720] - rev_rates[671]);

  //rxn 721
  //sp 4
  sp_rates[3] += (fwd_rates[721] - rev_rates[672]);
  //sp 106
  sp_rates[105] = (fwd_rates[721] - rev_rates[672]);
  //sp 36
  sp_rates[35] -= (fwd_rates[721] - rev_rates[672]);
  //sp 37
  sp_rates[36] -= (fwd_rates[721] - rev_rates[672]);

  //rxn 722
  //sp 106
  sp_rates[105] -= (fwd_rates[722] - rev_rates[673]);
  //sp 11
  sp_rates[10] += (fwd_rates[722] - rev_rates[673]);
  //sp 80
  sp_rates[79] += (fwd_rates[722] - rev_rates[673]);
  //sp 8
  sp_rates[7] -= (fwd_rates[722] - rev_rates[673]);

  //rxn 723
  //sp 107
  sp_rates[106] += (fwd_rates[723] - rev_rates[674]) * pres_mod[42];
  //sp 4
  sp_rates[3] -= (fwd_rates[723] - rev_rates[674]) * pres_mod[42];
  //sp 109
  sp_rates[108] -= (fwd_rates[723] - rev_rates[674]) * pres_mod[42];

  //rxn 724
  //sp 107
  sp_rates[106] += (fwd_rates[724] - rev_rates[675]);
  //sp 108
  sp_rates[107] -= (fwd_rates[724] - rev_rates[675]);

  //rxn 725
  //sp 4
  sp_rates[3] += (fwd_rates[725] - rev_rates[676]);
  //sp 108
  sp_rates[107] -= (fwd_rates[725] - rev_rates[676]);
  //sp 109
  sp_rates[108] += (fwd_rates[725] - rev_rates[676]);

  //rxn 726
  //sp 108
  sp_rates[107] += (fwd_rates[726] - rev_rates[677]);
  //sp 107
  sp_rates[106] -= (fwd_rates[726] - rev_rates[677]);

  //rxn 727
  //sp 3
  sp_rates[2] += fwd_rates[727];
  //sp 107
  sp_rates[106] -= fwd_rates[727];
  //sp 4
  sp_rates[3] -= fwd_rates[727];
  //sp 109
  sp_rates[108] += fwd_rates[727];

  //rxn 728
  //sp 3
  sp_rates[2] -= fwd_rates[728];
  //sp 4
  sp_rates[3] += fwd_rates[728];
  //sp 109
  sp_rates[108] -= fwd_rates[728];
  //sp 107
  sp_rates[106] += fwd_rates[728];

  //rxn 729
  //sp 107
  sp_rates[106] -= (fwd_rates[729] - rev_rates[678]);
  //sp 109
  sp_rates[108] += (fwd_rates[729] - rev_rates[678]);
  //sp 7
  sp_rates[6] += (fwd_rates[729] - rev_rates[678]);
  //sp 8
  sp_rates[7] -= (fwd_rates[729] - rev_rates[678]);

  //rxn 730
  //sp 109
  sp_rates[108] += (fwd_rates[730] - rev_rates[679]);
  //sp 10
  sp_rates[9] += (fwd_rates[730] - rev_rates[679]);
  //sp 107
  sp_rates[106] -= (fwd_rates[730] - rev_rates[679]);
  //sp 5
  sp_rates[4] -= (fwd_rates[730] - rev_rates[679]);

  //rxn 731
  //sp 107
  sp_rates[106] -= (fwd_rates[731] - rev_rates[680]);
  //sp 109
  sp_rates[108] += (fwd_rates[731] - rev_rates[680]);
  //sp 6
  sp_rates[5] -= (fwd_rates[731] - rev_rates[680]);
  //sp 8
  sp_rates[7] += (fwd_rates[731] - rev_rates[680]);

  //rxn 732
  //sp 9
  sp_rates[8] += (fwd_rates[732] - rev_rates[681]);
  //sp 10
  sp_rates[9] -= (fwd_rates[732] - rev_rates[681]);
  //sp 107
  sp_rates[106] -= (fwd_rates[732] - rev_rates[681]);
  //sp 109
  sp_rates[108] += (fwd_rates[732] - rev_rates[681]);

  //rxn 733
  //sp 25
  sp_rates[24] += (fwd_rates[733] - rev_rates[682]);
  //sp 26
  sp_rates[25] -= (fwd_rates[733] - rev_rates[682]);
  //sp 107
  sp_rates[106] -= (fwd_rates[733] - rev_rates[682]);
  //sp 109
  sp_rates[108] += (fwd_rates[733] - rev_rates[682]);

  //rxn 734
  //sp 108
  sp_rates[107] -= fwd_rates[734];
  //sp 3
  sp_rates[2] += fwd_rates[734];
  //sp 4
  sp_rates[3] -= fwd_rates[734];
  //sp 109
  sp_rates[108] += fwd_rates[734];

  //rxn 735
  //sp 108
  sp_rates[107] -= fwd_rates[735];
  //sp 109
  sp_rates[108] += fwd_rates[735];
  //sp 7
  sp_rates[6] += fwd_rates[735];
  //sp 8
  sp_rates[7] -= fwd_rates[735];

  //rxn 736
  //sp 108
  sp_rates[107] -= fwd_rates[736];
  //sp 109
  sp_rates[108] += fwd_rates[736];
  //sp 6
  sp_rates[5] -= fwd_rates[736];
  //sp 8
  sp_rates[7] += fwd_rates[736];

  //rxn 737
  //sp 109
  sp_rates[108] += fwd_rates[737];
  //sp 10
  sp_rates[9] += fwd_rates[737];
  //sp 108
  sp_rates[107] -= fwd_rates[737];
  //sp 5
  sp_rates[4] -= fwd_rates[737];

  //rxn 738
  //sp 9
  sp_rates[8] += fwd_rates[738];
  //sp 10
  sp_rates[9] -= fwd_rates[738];
  //sp 108
  sp_rates[107] -= fwd_rates[738];
  //sp 109
  sp_rates[108] += fwd_rates[738];

  //rxn 739
  //sp 25
  sp_rates[24] += fwd_rates[739];
  //sp 26
  sp_rates[25] -= fwd_rates[739];
  //sp 108
  sp_rates[107] -= fwd_rates[739];
  //sp 109
  sp_rates[108] += fwd_rates[739];

  //rxn 740
  //sp 108
  sp_rates[107] -= fwd_rates[740];
  //sp 109
  sp_rates[108] += fwd_rates[740];
  //sp 14
  sp_rates[13] += fwd_rates[740];
  //sp 15
  sp_rates[14] -= fwd_rates[740];

  //rxn 741
  //sp 26
  sp_rates[25] -= (fwd_rates[741] - rev_rates[683]);
  //sp 11
  sp_rates[10] += (fwd_rates[741] - rev_rates[683]);
  //sp 109
  sp_rates[108] -= (fwd_rates[741] - rev_rates[683]);
  //sp 107
  sp_rates[106] += (fwd_rates[741] - rev_rates[683]);

  //rxn 742
  //sp 36
  sp_rates[35] += fwd_rates[742];
  //sp 4
  sp_rates[3] += fwd_rates[742];
  //sp 109
  sp_rates[108] -= fwd_rates[742];
  //sp 106
  sp_rates[105] += fwd_rates[742];

  //rxn 743
  //sp 4
  sp_rates[3] += fwd_rates[743];
  //sp 102
  sp_rates[101] -= fwd_rates[743];
  //sp 103
  sp_rates[102] = fwd_rates[743];

  //rxn 744
  //sp 4
  sp_rates[3] += (fwd_rates[744] - rev_rates[684]);
  //sp 102
  sp_rates[101] -= (fwd_rates[744] - rev_rates[684]);
  //sp 103
  sp_rates[102] += (fwd_rates[744] - rev_rates[684]);

  //rxn 745
  //sp 10
  sp_rates[9] -= (fwd_rates[745] - rev_rates[685]);
  //sp 5
  sp_rates[4] += (fwd_rates[745] - rev_rates[685]);
  //sp 102
  sp_rates[101] += (fwd_rates[745] - rev_rates[685]);
  //sp 103
  sp_rates[102] -= (fwd_rates[745] - rev_rates[685]);

  //rxn 746
  //sp 102
  sp_rates[101] -= (fwd_rates[746] - rev_rates[686]);
  //sp 6
  sp_rates[5] -= (fwd_rates[746] - rev_rates[686]);
  //sp 103
  sp_rates[102] += (fwd_rates[746] - rev_rates[686]);
  //sp 8
  sp_rates[7] += (fwd_rates[746] - rev_rates[686]);

  //rxn 747
  //sp 103
  sp_rates[102] += (fwd_rates[747] - rev_rates[687]);
  //sp 102
  sp_rates[101] -= (fwd_rates[747] - rev_rates[687]);
  //sp 7
  sp_rates[6] += (fwd_rates[747] - rev_rates[687]);
  //sp 8
  sp_rates[7] -= (fwd_rates[747] - rev_rates[687]);

  //rxn 748
  //sp 9
  sp_rates[8] -= (fwd_rates[748] - rev_rates[688]);
  //sp 10
  sp_rates[9] += (fwd_rates[748] - rev_rates[688]);
  //sp 102
  sp_rates[101] += (fwd_rates[748] - rev_rates[688]);
  //sp 103
  sp_rates[102] -= (fwd_rates[748] - rev_rates[688]);

  //rxn 749
  //sp 34
  sp_rates[33] += (fwd_rates[749] - rev_rates[689]);
  //sp 35
  sp_rates[34] -= (fwd_rates[749] - rev_rates[689]);
  //sp 102
  sp_rates[101] -= (fwd_rates[749] - rev_rates[689]);
  //sp 103
  sp_rates[102] += (fwd_rates[749] - rev_rates[689]);

  //rxn 750
  //sp 103
  sp_rates[102] += (fwd_rates[750] - rev_rates[690]);
  //sp 102
  sp_rates[101] -= (fwd_rates[750] - rev_rates[690]);
  //sp 79
  sp_rates[78] += (fwd_rates[750] - rev_rates[690]);
  //sp 80
  sp_rates[79] -= (fwd_rates[750] - rev_rates[690]);

  //rxn 751
  //sp 103
  sp_rates[102] += (fwd_rates[751] - rev_rates[691]);
  //sp 62
  sp_rates[61] += (fwd_rates[751] - rev_rates[691]);
  //sp 102
  sp_rates[101] -= (fwd_rates[751] - rev_rates[691]);
  //sp 63
  sp_rates[62] -= (fwd_rates[751] - rev_rates[691]);

  //rxn 752
  //sp 35
  sp_rates[34] -= (fwd_rates[752] - rev_rates[692]);
  //sp 36
  sp_rates[35] -= (fwd_rates[752] - rev_rates[692]);
  //sp 103
  sp_rates[102] += (fwd_rates[752] - rev_rates[692]);

  //rxn 753
  //sp 35
  sp_rates[34] -= 2.0 * (fwd_rates[753] - rev_rates[693]);
  //sp 4
  sp_rates[3] += (fwd_rates[753] - rev_rates[693]);
  //sp 103
  sp_rates[102] += (fwd_rates[753] - rev_rates[693]);

  //rxn 754
  //sp 4
  sp_rates[3] -= (fwd_rates[754] - rev_rates[694]);
  //sp 103
  sp_rates[102] += (fwd_rates[754] - rev_rates[694]);
  //sp 104
  sp_rates[103] -= (fwd_rates[754] - rev_rates[694]);

  //rxn 755
  //sp 4
  sp_rates[3] -= (fwd_rates[755] - rev_rates[695]);
  //sp 103
  sp_rates[102] += (fwd_rates[755] - rev_rates[695]);
  //sp 104
  sp_rates[103] -= (fwd_rates[755] - rev_rates[695]);

  //rxn 756
  //sp 3
  sp_rates[2] += (fwd_rates[756] - rev_rates[696]);
  //sp 4
  sp_rates[3] -= (fwd_rates[756] - rev_rates[696]);
  //sp 103
  sp_rates[102] -= (fwd_rates[756] - rev_rates[696]);
  //sp 104
  sp_rates[103] += (fwd_rates[756] - rev_rates[696]);

  //rxn 757
  //sp 7
  sp_rates[6] += (fwd_rates[757] - rev_rates[697]);
  //sp 104
  sp_rates[103] += (fwd_rates[757] - rev_rates[697]);
  //sp 103
  sp_rates[102] -= (fwd_rates[757] - rev_rates[697]);
  //sp 8
  sp_rates[7] -= (fwd_rates[757] - rev_rates[697]);

  //rxn 758
  //sp 26
  sp_rates[25] -= (fwd_rates[758] - rev_rates[698]);
  //sp 11
  sp_rates[10] += (fwd_rates[758] - rev_rates[698]);
  //sp 102
  sp_rates[101] += (fwd_rates[758] - rev_rates[698]);
  //sp 103
  sp_rates[102] -= (fwd_rates[758] - rev_rates[698]);

  //rxn 759
  //sp 35
  sp_rates[34] += fwd_rates[759];
  //sp 103
  sp_rates[102] -= fwd_rates[759];
  //sp 8
  sp_rates[7] += fwd_rates[759];
  //sp 10
  sp_rates[9] -= fwd_rates[759];
  //sp 48
  sp_rates[47] += fwd_rates[759];

  //rxn 760
  //sp 26
  sp_rates[25] += (fwd_rates[760] - rev_rates[699]);
  //sp 5
  sp_rates[4] -= (fwd_rates[760] - rev_rates[699]);
  //sp 103
  sp_rates[102] -= (fwd_rates[760] - rev_rates[699]);
  //sp 85
  sp_rates[84] += (fwd_rates[760] - rev_rates[699]);

  //rxn 761
  //sp 4
  sp_rates[3] += fwd_rates[761];
  //sp 5
  sp_rates[4] -= fwd_rates[761];
  //sp 103
  sp_rates[102] -= fwd_rates[761];
  //sp 11
  sp_rates[10] += fwd_rates[761];
  //sp 85
  sp_rates[84] += fwd_rates[761];

  //rxn 762
  //sp 4
  sp_rates[3] += (fwd_rates[762] - rev_rates[700]);
  //sp 108
  sp_rates[107] += (fwd_rates[762] - rev_rates[700]);
  //sp 36
  sp_rates[35] -= (fwd_rates[762] - rev_rates[700]);
  //sp 103
  sp_rates[102] -= (fwd_rates[762] - rev_rates[700]);

  //rxn 763
  //sp 4
  sp_rates[3] += (fwd_rates[763] - rev_rates[701]);
  //sp 107
  sp_rates[106] += (fwd_rates[763] - rev_rates[701]);
  //sp 36
  sp_rates[35] -= (fwd_rates[763] - rev_rates[701]);
  //sp 103
  sp_rates[102] -= (fwd_rates[763] - rev_rates[701]);

  //rxn 764
  //sp 107
  sp_rates[106] += (fwd_rates[764] - rev_rates[702]);
  //sp 36
  sp_rates[35] -= (fwd_rates[764] - rev_rates[702]);
  //sp 104
  sp_rates[103] -= (fwd_rates[764] - rev_rates[702]);

  //rxn 765
  //sp 25
  sp_rates[24] += (fwd_rates[765] - rev_rates[703]);
  //sp 102
  sp_rates[101] -= (fwd_rates[765] - rev_rates[703]);
  //sp 6
  sp_rates[5] -= (fwd_rates[765] - rev_rates[703]);
  //sp 79
  sp_rates[78] += (fwd_rates[765] - rev_rates[703]);

  //rxn 766
  //sp 102
  sp_rates[101] -= (fwd_rates[766] - rev_rates[704]);
  //sp 35
  sp_rates[34] += (fwd_rates[766] - rev_rates[704]);
  //sp 6
  sp_rates[5] -= (fwd_rates[766] - rev_rates[704]);
  //sp 47
  sp_rates[46] += (fwd_rates[766] - rev_rates[704]);

  //rxn 767
  //sp 35
  sp_rates[34] -= (fwd_rates[767] - rev_rates[705]);
  //sp 3
  sp_rates[2] += (fwd_rates[767] - rev_rates[705]);
  //sp 103
  sp_rates[102] -= (fwd_rates[767] - rev_rates[705]);
  //sp 107
  sp_rates[106] += (fwd_rates[767] - rev_rates[705]);

  //rxn 768
  //sp 37
  sp_rates[36] -= (fwd_rates[768] - rev_rates[706]);
  //sp 103
  sp_rates[102] -= (fwd_rates[768] - rev_rates[706]);
  //sp 80
  sp_rates[79] += 2.0 * (fwd_rates[768] - rev_rates[706]);

  //rxn 769
  //sp 105
  sp_rates[104] = -(fwd_rates[769] - rev_rates[707]);
  //sp 4
  sp_rates[3] -= (fwd_rates[769] - rev_rates[707]);
  //sp 104
  sp_rates[103] += (fwd_rates[769] - rev_rates[707]);

  //rxn 770
  //sp 105
  sp_rates[104] += (fwd_rates[770] - rev_rates[708]);
  //sp 3
  sp_rates[2] += (fwd_rates[770] - rev_rates[708]);
  //sp 4
  sp_rates[3] -= (fwd_rates[770] - rev_rates[708]);
  //sp 104
  sp_rates[103] -= (fwd_rates[770] - rev_rates[708]);

  //rxn 771
  //sp 105
  sp_rates[104] += (fwd_rates[771] - rev_rates[709]);
  //sp 104
  sp_rates[103] -= (fwd_rates[771] - rev_rates[709]);
  //sp 7
  sp_rates[6] += (fwd_rates[771] - rev_rates[709]);
  //sp 8
  sp_rates[7] -= (fwd_rates[771] - rev_rates[709]);

  //rxn 772
  //sp 105
  sp_rates[104] += (fwd_rates[772] - rev_rates[710]);
  //sp 14
  sp_rates[13] += (fwd_rates[772] - rev_rates[710]);
  //sp 15
  sp_rates[14] -= (fwd_rates[772] - rev_rates[710]);
  //sp 104
  sp_rates[103] -= (fwd_rates[772] - rev_rates[710]);

  //rxn 773
  //sp 105
  sp_rates[104] += (fwd_rates[773] - rev_rates[711]) * pres_mod[43];
  //sp 36
  sp_rates[35] -= (fwd_rates[773] - rev_rates[711]) * pres_mod[43];
  //sp 37
  sp_rates[36] -= (fwd_rates[773] - rev_rates[711]) * pres_mod[43];

  //rxn 774
  //sp 105
  sp_rates[104] += (fwd_rates[774] - rev_rates[712]);
  //sp 106
  sp_rates[105] -= (fwd_rates[774] - rev_rates[712]);
  //sp 4
  sp_rates[3] -= (fwd_rates[774] - rev_rates[712]);

  //rxn 775
  //sp 105
  sp_rates[104] -= (fwd_rates[775] - rev_rates[713]);
  //sp 106
  sp_rates[105] += (fwd_rates[775] - rev_rates[713]);
  //sp 3
  sp_rates[2] += (fwd_rates[775] - rev_rates[713]);
  //sp 4
  sp_rates[3] -= (fwd_rates[775] - rev_rates[713]);

  //rxn 776
  //sp 105
  sp_rates[104] -= (fwd_rates[776] - rev_rates[714]);
  //sp 106
  sp_rates[105] += (fwd_rates[776] - rev_rates[714]);
  //sp 7
  sp_rates[6] += (fwd_rates[776] - rev_rates[714]);
  //sp 8
  sp_rates[7] -= (fwd_rates[776] - rev_rates[714]);

  //rxn 777
  //sp 105
  sp_rates[104] -= (fwd_rates[777] - rev_rates[715]);
  //sp 35
  sp_rates[34] -= (fwd_rates[777] - rev_rates[715]);
  //sp 80
  sp_rates[79] += 2.0 * (fwd_rates[777] - rev_rates[715]);

  //rxn 778
  //sp 105
  sp_rates[104] += (fwd_rates[778] - rev_rates[716]);
  //sp 19
  sp_rates[18] -= (fwd_rates[778] - rev_rates[716]);
  //sp 4
  sp_rates[3] += (fwd_rates[778] - rev_rates[716]);
  //sp 80
  sp_rates[79] -= (fwd_rates[778] - rev_rates[716]);

  //rxn 779
  //sp 105
  sp_rates[104] -= (fwd_rates[779] - rev_rates[717]);
  //sp 49
  sp_rates[48] += (fwd_rates[779] - rev_rates[717]);
  //sp 5
  sp_rates[4] -= (fwd_rates[779] - rev_rates[717]);
  //sp 48
  sp_rates[47] += (fwd_rates[779] - rev_rates[717]);

  //rxn 780
  //sp 105
  sp_rates[104] -= (fwd_rates[780] - rev_rates[718]);
  //sp 37
  sp_rates[36] += (fwd_rates[780] - rev_rates[718]);
  //sp 79
  sp_rates[78] += (fwd_rates[780] - rev_rates[718]);
  //sp 16
  sp_rates[15] -= (fwd_rates[780] - rev_rates[718]);

  //rxn 781
  //sp 3
  sp_rates[2] += (fwd_rates[781] - rev_rates[719]);
  //sp 45
  sp_rates[44] -= (fwd_rates[781] - rev_rates[719]);
  //sp 48
  sp_rates[47] += (fwd_rates[781] - rev_rates[719]);

  //rxn 782
  //sp 25
  sp_rates[24] += fwd_rates[782];
  //sp 5
  sp_rates[4] -= fwd_rates[782];
  //sp 79
  sp_rates[78] -= fwd_rates[782];
  //sp 48
  sp_rates[47] += fwd_rates[782];

  //rxn 783
  //sp 36
  sp_rates[35] -= fwd_rates[783];
  //sp 6
  sp_rates[5] -= fwd_rates[783];
  //sp 48
  sp_rates[47] += fwd_rates[783];

  //rxn 784
  //sp 26
  sp_rates[25] += 2.0 * fwd_rates[784];
  //sp 6
  sp_rates[5] -= fwd_rates[784];
  //sp 48
  sp_rates[47] -= fwd_rates[784];

  //rxn 785
  //sp 25
  sp_rates[24] += fwd_rates[785];
  //sp 11
  sp_rates[10] += fwd_rates[785];
  //sp 6
  sp_rates[5] -= fwd_rates[785];
  //sp 48
  sp_rates[47] -= fwd_rates[785];

  //rxn 786
  //sp 12
  sp_rates[11] += fwd_rates[786];
  //sp 48
  sp_rates[47] -= fwd_rates[786];
  //sp 15
  sp_rates[14] += fwd_rates[786];
  //sp 8
  sp_rates[7] -= fwd_rates[786];

  //rxn 787
  //sp 25
  sp_rates[24] += fwd_rates[787];
  //sp 12
  sp_rates[11] += fwd_rates[787];
  //sp 5
  sp_rates[4] -= fwd_rates[787];
  //sp 48
  sp_rates[47] -= fwd_rates[787];

  //rxn 788
  //sp 5
  sp_rates[4] -= fwd_rates[788];
  //sp 8
  sp_rates[7] += fwd_rates[788];
  //sp 11
  sp_rates[10] += fwd_rates[788];
  //sp 48
  sp_rates[47] -= fwd_rates[788];
  //sp 26
  sp_rates[25] += fwd_rates[788];

  //rxn 789
  //sp 105
  sp_rates[104] -= (fwd_rates[789] - rev_rates[720]);
  //sp 37
  sp_rates[36] += (fwd_rates[789] - rev_rates[720]);
  //sp 6
  sp_rates[5] -= (fwd_rates[789] - rev_rates[720]);
  //sp 48
  sp_rates[47] += (fwd_rates[789] - rev_rates[720]);

  //rxn 790
  //sp 36
  sp_rates[35] += fwd_rates[790];
  //sp 5
  sp_rates[4] -= fwd_rates[790];
  //sp 103
  sp_rates[102] -= fwd_rates[790];
  //sp 8
  sp_rates[7] += fwd_rates[790];
  //sp 48
  sp_rates[47] += fwd_rates[790];

  //rxn 791
  //sp 34
  sp_rates[33] += fwd_rates[791];
  //sp 48
  sp_rates[47] += fwd_rates[791];
  //sp 103
  sp_rates[102] -= fwd_rates[791];
  //sp 8
  sp_rates[7] -= fwd_rates[791];

  //rxn 792
  //sp 48
  sp_rates[47] += fwd_rates[792];
  //sp 46
  sp_rates[45] -= fwd_rates[792];
  //sp 7
  sp_rates[6] += fwd_rates[792];
  //sp 8
  sp_rates[7] -= fwd_rates[792];

  //rxn 793
  //sp 7
  sp_rates[6] += fwd_rates[793];
  //sp 48
  sp_rates[47] += fwd_rates[793];
  //sp 47
  sp_rates[46] -= fwd_rates[793];
  //sp 8
  sp_rates[7] -= fwd_rates[793];

  //rxn 794
  //sp 8
  sp_rates[7] += fwd_rates[794];
  //sp 10
  sp_rates[9] -= fwd_rates[794];
  //sp 11
  sp_rates[10] += fwd_rates[794];
  //sp 48
  sp_rates[47] -= fwd_rates[794];
  //sp 25
  sp_rates[24] += fwd_rates[794];

  //rxn 795
  //sp 10
  sp_rates[9] += fwd_rates[795];
  //sp 5
  sp_rates[4] -= fwd_rates[795];
  //sp 47
  sp_rates[46] -= fwd_rates[795];
  //sp 48
  sp_rates[47] += fwd_rates[795];

  //rxn 796
  //sp 49
  sp_rates[48] += (fwd_rates[796] - rev_rates[721]);
  //sp 46
  sp_rates[45] += (fwd_rates[796] - rev_rates[721]);
  //sp 48
  sp_rates[47] -= 2.0 * (fwd_rates[796] - rev_rates[721]);

  //rxn 797
  //sp 12
  sp_rates[11] += (fwd_rates[797] - rev_rates[722]);
  //sp 79
  sp_rates[78] += (fwd_rates[797] - rev_rates[722]);
  //sp 48
  sp_rates[47] -= 2.0 * (fwd_rates[797] - rev_rates[722]);

  //rxn 798
  //sp 34
  sp_rates[33] += fwd_rates[798];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[798];
  //sp 48
  sp_rates[47] -= 2.0 * fwd_rates[798];

  //rxn 799
  //sp 49
  sp_rates[48] += fwd_rates[799];
  //sp 4
  sp_rates[3] += fwd_rates[799];
  //sp 48
  sp_rates[47] -= fwd_rates[799];

  //rxn 800
  //sp 34
  sp_rates[33] += (fwd_rates[800] - rev_rates[723]);
  //sp 11
  sp_rates[10] += (fwd_rates[800] - rev_rates[723]);
  //sp 48
  sp_rates[47] -= (fwd_rates[800] - rev_rates[723]);
  //sp 16
  sp_rates[15] -= (fwd_rates[800] - rev_rates[723]);

  //rxn 801
  //sp 49
  sp_rates[48] += (fwd_rates[801] - rev_rates[724]);
  //sp 48
  sp_rates[47] -= (fwd_rates[801] - rev_rates[724]);
  //sp 15
  sp_rates[14] += (fwd_rates[801] - rev_rates[724]);
  //sp 16
  sp_rates[15] -= (fwd_rates[801] - rev_rates[724]);

  //rxn 802
  //sp 79
  sp_rates[78] += (fwd_rates[802] - rev_rates[725]);
  //sp 8
  sp_rates[7] += (fwd_rates[802] - rev_rates[725]);
  //sp 15
  sp_rates[14] -= (fwd_rates[802] - rev_rates[725]);
  //sp 48
  sp_rates[47] -= (fwd_rates[802] - rev_rates[725]);

  //rxn 803
  //sp 8
  sp_rates[7] += (fwd_rates[803] - rev_rates[726]);
  //sp 78
  sp_rates[77] += (fwd_rates[803] - rev_rates[726]);
  //sp 15
  sp_rates[14] -= (fwd_rates[803] - rev_rates[726]);
  //sp 48
  sp_rates[47] -= (fwd_rates[803] - rev_rates[726]);

  //rxn 804
  //sp 4
  sp_rates[3] += (fwd_rates[804] - rev_rates[727]);
  //sp 8
  sp_rates[7] -= (fwd_rates[804] - rev_rates[727]);
  //sp 11
  sp_rates[10] += (fwd_rates[804] - rev_rates[727]);
  //sp 49
  sp_rates[48] -= (fwd_rates[804] - rev_rates[727]);
  //sp 26
  sp_rates[25] += (fwd_rates[804] - rev_rates[727]);

  //rxn 805
  //sp 49
  sp_rates[48] -= (fwd_rates[805] - rev_rates[728]);
  //sp 26
  sp_rates[25] += (fwd_rates[805] - rev_rates[728]);
  //sp 12
  sp_rates[11] += (fwd_rates[805] - rev_rates[728]);
  //sp 5
  sp_rates[4] -= (fwd_rates[805] - rev_rates[728]);

  //rxn 806
  //sp 49
  sp_rates[48] -= (fwd_rates[806] - rev_rates[729]);
  //sp 11
  sp_rates[10] += (fwd_rates[806] - rev_rates[729]);
  //sp 35
  sp_rates[34] += (fwd_rates[806] - rev_rates[729]);
  //sp 16
  sp_rates[15] -= (fwd_rates[806] - rev_rates[729]);

  //rxn 807
  //sp 49
  sp_rates[48] -= 2.0 * fwd_rates[807];
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[807];
  //sp 36
  sp_rates[35] += fwd_rates[807];

  //rxn 808
  //sp 49
  sp_rates[48] -= (fwd_rates[808] - rev_rates[730]);
  //sp 35
  sp_rates[34] -= (fwd_rates[808] - rev_rates[730]);
  //sp 11
  sp_rates[10] += (fwd_rates[808] - rev_rates[730]);
  //sp 78
  sp_rates[77] += (fwd_rates[808] - rev_rates[730]);

  //rxn 809
  //sp 49
  sp_rates[48] -= (fwd_rates[809] - rev_rates[731]);
  //sp 11
  sp_rates[10] += (fwd_rates[809] - rev_rates[731]);
  //sp 102
  sp_rates[101] += (fwd_rates[809] - rev_rates[731]);
  //sp 63
  sp_rates[62] -= (fwd_rates[809] - rev_rates[731]);

  //rxn 810
  //sp 49
  sp_rates[48] += (fwd_rates[810] - rev_rates[732]);
  //sp 36
  sp_rates[35] -= (fwd_rates[810] - rev_rates[732]);
  //sp 5
  sp_rates[4] -= (fwd_rates[810] - rev_rates[732]);
  //sp 8
  sp_rates[7] += (fwd_rates[810] - rev_rates[732]);

  //rxn 811
  //sp 49
  sp_rates[48] += (fwd_rates[811] - rev_rates[733]);
  //sp 106
  sp_rates[105] -= (fwd_rates[811] - rev_rates[733]);
  //sp 36
  sp_rates[35] += (fwd_rates[811] - rev_rates[733]);
  //sp 8
  sp_rates[7] -= (fwd_rates[811] - rev_rates[733]);

  //rxn 812
  //sp 49
  sp_rates[48] += (fwd_rates[812] - rev_rates[734]);
  //sp 5
  sp_rates[4] -= (fwd_rates[812] - rev_rates[734]);
  //sp 6
  sp_rates[5] += (fwd_rates[812] - rev_rates[734]);
  //sp 37
  sp_rates[36] -= (fwd_rates[812] - rev_rates[734]);

  //rxn 813
  //sp 81
  sp_rates[80] -= (fwd_rates[813] - rev_rates[735]);
  //sp 26
  sp_rates[25] += (fwd_rates[813] - rev_rates[735]);
  //sp 49
  sp_rates[48] += (fwd_rates[813] - rev_rates[735]);
  //sp 5
  sp_rates[4] -= (fwd_rates[813] - rev_rates[735]);

  //rxn 814
  //sp 85
  sp_rates[84] += (fwd_rates[814] - rev_rates[736]);
  //sp 80
  sp_rates[79] -= (fwd_rates[814] - rev_rates[736]);
  //sp 8
  sp_rates[7] -= (fwd_rates[814] - rev_rates[736]);

  //rxn 815
  //sp 4
  sp_rates[3] += fwd_rates[815];
  //sp 85
  sp_rates[84] += fwd_rates[815];
  //sp 78
  sp_rates[77] -= fwd_rates[815];
  //sp 8
  sp_rates[7] -= fwd_rates[815];

  //rxn 816
  //sp 26
  sp_rates[25] += fwd_rates[816];
  //sp 85
  sp_rates[84] -= fwd_rates[816];
  //sp 45
  sp_rates[44] += fwd_rates[816];
  //sp 8
  sp_rates[7] -= fwd_rates[816];

  //rxn 817
  //sp 34
  sp_rates[33] += fwd_rates[817];
  //sp 4
  sp_rates[3] += fwd_rates[817];
  //sp 8
  sp_rates[7] -= fwd_rates[817];
  //sp 12
  sp_rates[11] += fwd_rates[817];
  //sp 85
  sp_rates[84] -= fwd_rates[817];

  //rxn 818
  //sp 8
  sp_rates[7] += fwd_rates[818];
  //sp 10
  sp_rates[9] -= fwd_rates[818];
  //sp 48
  sp_rates[47] += fwd_rates[818];
  //sp 85
  sp_rates[84] -= fwd_rates[818];
  //sp 25
  sp_rates[24] += fwd_rates[818];

  //rxn 819
  //sp 34
  sp_rates[33] -= fwd_rates[819];
  //sp 4
  sp_rates[3] += 0.2 * fwd_rates[819];
  //sp 11
  sp_rates[10] += 0.8 * fwd_rates[819];
  //sp 85
  sp_rates[84] += 0.2 * fwd_rates[819];
  //sp 26
  sp_rates[25] -= fwd_rates[819];
  //sp 31
  sp_rates[30] += 0.8 * fwd_rates[819];

  //rxn 820
  //sp 26
  sp_rates[25] -= fwd_rates[820];
  //sp 85
  sp_rates[84] += fwd_rates[820];
  //sp 62
  sp_rates[61] -= fwd_rates[820];
  //sp 15
  sp_rates[14] += fwd_rates[820];

  //rxn 821
  //sp 25
  sp_rates[24] += fwd_rates[821];
  //sp 5
  sp_rates[4] -= fwd_rates[821];
  //sp 102
  sp_rates[101] -= fwd_rates[821];
  //sp 85
  sp_rates[84] += fwd_rates[821];

  //rxn 822
  //sp 85
  sp_rates[84] += fwd_rates[822];
  //sp 5
  sp_rates[4] -= fwd_rates[822];
  //sp 110
  sp_rates[109] = -fwd_rates[822];
  //sp 48
  sp_rates[47] += fwd_rates[822];

  //rxn 823
  //sp 26
  sp_rates[25] += fwd_rates[823];
  //sp 4
  sp_rates[3] -= fwd_rates[823];
  //sp 85
  sp_rates[84] -= fwd_rates[823];
  //sp 34
  sp_rates[33] += fwd_rates[823];

  //rxn 824
  //sp 47
  sp_rates[46] += fwd_rates[824];
  //sp 67
  sp_rates[66] = -fwd_rates[824];
  //sp 15
  sp_rates[14] += fwd_rates[824];

  //rxn 825
  //sp 25
  sp_rates[24] += fwd_rates[825];
  //sp 34
  sp_rates[33] += fwd_rates[825];
  //sp 67
  sp_rates[66] -= fwd_rates[825];

  //rxn 826
  //sp 67
  sp_rates[66] += fwd_rates[826];
  //sp 62
  sp_rates[61] -= fwd_rates[826];
  //sp 6
  sp_rates[5] -= fwd_rates[826];

  //rxn 827
  //sp 34
  sp_rates[33] += 0.5 * fwd_rates[827];
  //sp 67
  sp_rates[66] += 0.5 * fwd_rates[827];
  //sp 8
  sp_rates[7] -= fwd_rates[827];
  //sp 25
  sp_rates[24] += 0.5 * fwd_rates[827];
  //sp 63
  sp_rates[62] -= fwd_rates[827];

  //rxn 828
  //sp 34
  sp_rates[33] += 0.5 * fwd_rates[828];
  //sp 67
  sp_rates[66] += 0.5 * fwd_rates[828];
  //sp 8
  sp_rates[7] -= fwd_rates[828];
  //sp 25
  sp_rates[24] += 0.5 * fwd_rates[828];
  //sp 64
  sp_rates[63] -= fwd_rates[828];

  //rxn 829
  //sp 67
  sp_rates[66] += 0.2 * fwd_rates[829];
  //sp 10
  sp_rates[9] -= fwd_rates[829];
  //sp 86
  sp_rates[85] += 0.8 * fwd_rates[829];
  //sp 25
  sp_rates[24] += fwd_rates[829];
  //sp 95
  sp_rates[94] -= fwd_rates[829];

  //rxn 830
  //sp 10
  sp_rates[9] += (fwd_rates[830] - rev_rates[737]) * pres_mod[44];
  //sp 6
  sp_rates[5] -= (fwd_rates[830] - rev_rates[737]) * pres_mod[44];
  //sp 8
  sp_rates[7] -= (fwd_rates[830] - rev_rates[737]) * pres_mod[44];

  //rxn 831
  //sp 19
  sp_rates[18] += (fwd_rates[831] - rev_rates[738]);
  //sp 12
  sp_rates[11] += (fwd_rates[831] - rev_rates[738]);
  //sp 5
  sp_rates[4] -= (fwd_rates[831] - rev_rates[738]);
  //sp 37
  sp_rates[36] -= (fwd_rates[831] - rev_rates[738]);

  //rxn 832
  //sp 18
  sp_rates[17] += (fwd_rates[832] - rev_rates[739]);
  //sp 19
  sp_rates[18] -= (fwd_rates[832] - rev_rates[739]);
  //sp 7
  sp_rates[6] += (fwd_rates[832] - rev_rates[739]);
  //sp 8
  sp_rates[7] -= (fwd_rates[832] - rev_rates[739]);

  //rxn 833
  //sp 25
  sp_rates[24] += fwd_rates[833];
  //sp 5
  sp_rates[4] -= fwd_rates[833];
  //sp 31
  sp_rates[30] -= fwd_rates[833];
  //sp 23
  sp_rates[22] += fwd_rates[833];

  //rxn 834
  //sp 5
  sp_rates[4] -= fwd_rates[834];
  //sp 6
  sp_rates[5] += fwd_rates[834];
  //sp 15
  sp_rates[14] += fwd_rates[834];
  //sp 25
  sp_rates[24] += fwd_rates[834];
  //sp 31
  sp_rates[30] -= fwd_rates[834];

  //rxn 835
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[835];
  //sp 4
  sp_rates[3] += 0.35 * fwd_rates[835];
  //sp 6
  sp_rates[5] -= fwd_rates[835];
  //sp 8
  sp_rates[7] += 0.3 * fwd_rates[835];
  //sp 45
  sp_rates[44] += 0.35 * fwd_rates[835];
  //sp 15
  sp_rates[14] += 0.35 * fwd_rates[835];
  //sp 25
  sp_rates[24] += 0.35 * fwd_rates[835];
  //sp 31
  sp_rates[30] -= fwd_rates[835];

  //rxn 836
  //sp 8
  sp_rates[7] += fwd_rates[836];
  //sp 10
  sp_rates[9] -= fwd_rates[836];
  //sp 15
  sp_rates[14] += fwd_rates[836];
  //sp 25
  sp_rates[24] += fwd_rates[836];
  //sp 31
  sp_rates[30] -= fwd_rates[836];

  //rxn 837
  //sp 34
  sp_rates[33] -= fwd_rates[837];
  //sp 45
  sp_rates[44] += fwd_rates[837];
  //sp 6
  sp_rates[5] -= fwd_rates[837];

  //rxn 838
  //sp 41
  sp_rates[40] += (fwd_rates[838] - rev_rates[740]);
  //sp 4
  sp_rates[3] += (fwd_rates[838] - rev_rates[740]);
  //sp 38
  sp_rates[37] -= (fwd_rates[838] - rev_rates[740]);

  //rxn 839
  //sp 25
  sp_rates[24] += 2.0 * fwd_rates[839];
  //sp 5
  sp_rates[4] -= fwd_rates[839];
  //sp 8
  sp_rates[7] += fwd_rates[839];
  //sp 40
  sp_rates[39] -= fwd_rates[839];

  //rxn 840
  //sp 4
  sp_rates[3] -= (fwd_rates[840] - rev_rates[741]);
  //sp 45
  sp_rates[44] += (fwd_rates[840] - rev_rates[741]);
  //sp 86
  sp_rates[85] -= (fwd_rates[840] - rev_rates[741]);
  //sp 15
  sp_rates[14] += (fwd_rates[840] - rev_rates[741]);

  //rxn 841
  //sp 48
  sp_rates[47] += (fwd_rates[841] - rev_rates[742]);
  //sp 87
  sp_rates[86] -= (fwd_rates[841] - rev_rates[742]);
  //sp 15
  sp_rates[14] += (fwd_rates[841] - rev_rates[742]);

  //rxn 842
  //sp 5
  sp_rates[4] -= fwd_rates[842];
  //sp 8
  sp_rates[7] += fwd_rates[842];
  //sp 48
  sp_rates[47] += fwd_rates[842];
  //sp 87
  sp_rates[86] -= fwd_rates[842];
  //sp 25
  sp_rates[24] += fwd_rates[842];

  //rxn 843
  //sp 8
  sp_rates[7] += fwd_rates[843];
  //sp 10
  sp_rates[9] -= fwd_rates[843];
  //sp 46
  sp_rates[45] += fwd_rates[843];
  //sp 87
  sp_rates[86] -= fwd_rates[843];
  //sp 25
  sp_rates[24] += fwd_rates[843];

  //rxn 844
  //sp 26
  sp_rates[25] += (fwd_rates[844] - rev_rates[743]);
  //sp 53
  sp_rates[52] = -(fwd_rates[844] - rev_rates[743]);
  //sp 24
  sp_rates[23] += (fwd_rates[844] - rev_rates[743]);

  //rxn 845
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[845] - rev_rates[744]);
  //sp 53
  sp_rates[52] -= (fwd_rates[845] - rev_rates[744]);

  //rxn 846
  //sp 53
  sp_rates[52] -= (fwd_rates[846] - rev_rates[745]);
  //sp 7
  sp_rates[6] += (fwd_rates[846] - rev_rates[745]);
  //sp 48
  sp_rates[47] += (fwd_rates[846] - rev_rates[745]);

  //rxn 847
  //sp 85
  sp_rates[84] -= fwd_rates[847];
  //sp 54
  sp_rates[53] = fwd_rates[847];
  //sp 15
  sp_rates[14] += fwd_rates[847];
  //sp 8
  sp_rates[7] -= fwd_rates[847];

  //rxn 848
  //sp 26
  sp_rates[25] += 2.0 * (fwd_rates[848] - rev_rates[746]);
  //sp 54
  sp_rates[53] -= (fwd_rates[848] - rev_rates[746]);

  //rxn 849
  //sp 50
  sp_rates[49] = -(fwd_rates[849] - rev_rates[747]);
  //sp 11
  sp_rates[10] += (fwd_rates[849] - rev_rates[747]);
  //sp 7
  sp_rates[6] += (fwd_rates[849] - rev_rates[747]);

  //rxn 850
  //sp 50
  sp_rates[49] -= (fwd_rates[850] - rev_rates[748]);
  //sp 3
  sp_rates[2] += (fwd_rates[850] - rev_rates[748]);
  //sp 12
  sp_rates[11] += (fwd_rates[850] - rev_rates[748]);

  //rxn 851
  //sp 57
  sp_rates[56] = -(fwd_rates[851] - rev_rates[749]);
  //sp 11
  sp_rates[10] += (fwd_rates[851] - rev_rates[749]);
  //sp 22
  sp_rates[21] += (fwd_rates[851] - rev_rates[749]);

  //rxn 852
  //sp 57
  sp_rates[56] -= (fwd_rates[852] - rev_rates[750]);
  //sp 12
  sp_rates[11] += (fwd_rates[852] - rev_rates[750]);
  //sp 14
  sp_rates[13] += (fwd_rates[852] - rev_rates[750]);

  //rxn 853
  //sp 57
  sp_rates[56] -= (fwd_rates[853] - rev_rates[751]);
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[853] - rev_rates[751]);

  //rxn 854
  //sp 57
  sp_rates[56] -= fwd_rates[854];
  //sp 12
  sp_rates[11] += fwd_rates[854];
  //sp 4
  sp_rates[3] += fwd_rates[854];
  //sp 15
  sp_rates[14] += fwd_rates[854];

  //rxn 855
  //sp 57
  sp_rates[56] += (fwd_rates[855] - rev_rates[752]);
  //sp 58
  sp_rates[57] = -(fwd_rates[855] - rev_rates[752]);
  //sp 4
  sp_rates[3] -= (fwd_rates[855] - rev_rates[752]);

  //rxn 856
  //sp 57
  sp_rates[56] += (fwd_rates[856] - rev_rates[753]);
  //sp 26
  sp_rates[25] -= (fwd_rates[856] - rev_rates[753]);
  //sp 23
  sp_rates[22] -= (fwd_rates[856] - rev_rates[753]);

  //rxn 857
  //sp 3
  sp_rates[2] += fwd_rates[857];
  //sp 4
  sp_rates[3] -= fwd_rates[857];
  //sp 57
  sp_rates[56] -= fwd_rates[857];
  //sp 25
  sp_rates[24] += fwd_rates[857];
  //sp 26
  sp_rates[25] += fwd_rates[857];

  //rxn 858
  //sp 57
  sp_rates[56] -= fwd_rates[858];
  //sp 58
  sp_rates[57] += fwd_rates[858];
  //sp 3
  sp_rates[2] += fwd_rates[858];
  //sp 4
  sp_rates[3] -= fwd_rates[858];

  //rxn 859
  //sp 58
  sp_rates[57] -= (fwd_rates[859] - rev_rates[754]);
  //sp 12
  sp_rates[11] += (fwd_rates[859] - rev_rates[754]);
  //sp 15
  sp_rates[14] += (fwd_rates[859] - rev_rates[754]);

  //rxn 860
  //sp 58
  sp_rates[57] -= (fwd_rates[860] - rev_rates[755]);
  //sp 11
  sp_rates[10] += (fwd_rates[860] - rev_rates[755]);
  //sp 23
  sp_rates[22] += (fwd_rates[860] - rev_rates[755]);

  //rxn 861
  //sp 20
  sp_rates[19] += (fwd_rates[861] - rev_rates[756]);
  //sp 4
  sp_rates[3] -= (fwd_rates[861] - rev_rates[756]);
  //sp 21
  sp_rates[20] -= (fwd_rates[861] - rev_rates[756]);

  //rxn 862
  //sp 25
  sp_rates[24] += (fwd_rates[862] - rev_rates[757]);
  //sp 21
  sp_rates[20] -= (fwd_rates[862] - rev_rates[757]);
  //sp 8
  sp_rates[7] += (fwd_rates[862] - rev_rates[757]);

  //rxn 863
  //sp 10
  sp_rates[9] += fwd_rates[863];
  //sp 21
  sp_rates[20] -= fwd_rates[863];
  //sp 23
  sp_rates[22] += fwd_rates[863];
  //sp 8
  sp_rates[7] -= fwd_rates[863];

  //rxn 864
  //sp 5
  sp_rates[4] += fwd_rates[864];
  //sp 7
  sp_rates[6] += fwd_rates[864];
  //sp 10
  sp_rates[9] -= fwd_rates[864];
  //sp 21
  sp_rates[20] -= fwd_rates[864];
  //sp 25
  sp_rates[24] += fwd_rates[864];

  //rxn 865
  //sp 25
  sp_rates[24] += fwd_rates[865];
  //sp 20
  sp_rates[19] += fwd_rates[865];
  //sp 21
  sp_rates[20] -= fwd_rates[865];
  //sp 23
  sp_rates[22] -= fwd_rates[865];

  //rxn 866
  //sp 97
  sp_rates[96] += fwd_rates[866];
  //sp 21
  sp_rates[20] -= fwd_rates[866];
  //sp 5
  sp_rates[4] += fwd_rates[866];
  //sp 63
  sp_rates[62] -= fwd_rates[866];

  //rxn 867
  //sp 35
  sp_rates[34] += fwd_rates[867];
  //sp 21
  sp_rates[20] -= fwd_rates[867];
  //sp 23
  sp_rates[22] += fwd_rates[867];
  //sp 25
  sp_rates[24] += fwd_rates[867];
  //sp 63
  sp_rates[62] -= fwd_rates[867];

  //rxn 868
  //sp 8
  sp_rates[7] += fwd_rates[868];
  //sp 20
  sp_rates[19] -= fwd_rates[868];
  //sp 22
  sp_rates[21] += fwd_rates[868];
  //sp 23
  sp_rates[22] -= fwd_rates[868];
  //sp 25
  sp_rates[24] += fwd_rates[868];

  //rxn 869
  //sp 3
  sp_rates[2] += fwd_rates[869];
  //sp 8
  sp_rates[7] += fwd_rates[869];
  //sp 11
  sp_rates[10] += fwd_rates[869];
  //sp 21
  sp_rates[20] -= fwd_rates[869];

  //rxn 870
  //sp 11
  sp_rates[10] -= fwd_rates[870];
  //sp 12
  sp_rates[11] += fwd_rates[870];
  //sp 21
  sp_rates[20] -= fwd_rates[870];
  //sp 23
  sp_rates[22] += fwd_rates[870];

  //rxn 871
  //sp 10
  sp_rates[9] -= (fwd_rates[871] - rev_rates[758]);
  //sp 31
  sp_rates[30] -= (fwd_rates[871] - rev_rates[758]);
  //sp 32
  sp_rates[31] += (fwd_rates[871] - rev_rates[758]);

  //rxn 872
  //sp 42
  sp_rates[41] -= (fwd_rates[872] - rev_rates[759]);
  //sp 5
  sp_rates[4] -= (fwd_rates[872] - rev_rates[759]);
  //sp 55
  sp_rates[54] = (fwd_rates[872] - rev_rates[759]);

  //rxn 873
  //sp 56
  sp_rates[55] = fwd_rates[873];
  //sp 55
  sp_rates[54] -= fwd_rates[873];
  //sp 8
  sp_rates[7] += fwd_rates[873];

  //rxn 874
  //sp 25
  sp_rates[24] += fwd_rates[874];
  //sp 26
  sp_rates[25] += fwd_rates[874];
  //sp 8
  sp_rates[7] += fwd_rates[874];
  //sp 56
  sp_rates[55] -= fwd_rates[874];

  //rxn 875
  //sp 33
  sp_rates[32] -= fwd_rates[875];
  //sp 25
  sp_rates[24] += fwd_rates[875];
  //sp 23
  sp_rates[22] += fwd_rates[875];

  //rxn 876
  //sp 33
  sp_rates[32] -= fwd_rates[876];
  //sp 5
  sp_rates[4] += fwd_rates[876];
  //sp 15
  sp_rates[14] += fwd_rates[876];
  //sp 21
  sp_rates[20] -= fwd_rates[876];
  //sp 23
  sp_rates[22] += fwd_rates[876];
  //sp 25
  sp_rates[24] += fwd_rates[876];

  //rxn 877
  //sp 33
  sp_rates[32] -= 2.0 * fwd_rates[877];
  //sp 25
  sp_rates[24] += 2.0 * fwd_rates[877];
  //sp 5
  sp_rates[4] += fwd_rates[877];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[877];

  //rxn 878
  //sp 33
  sp_rates[32] -= fwd_rates[878];
  //sp 11
  sp_rates[10] -= fwd_rates[878];
  //sp 12
  sp_rates[11] += fwd_rates[878];
  //sp 15
  sp_rates[14] += fwd_rates[878];
  //sp 25
  sp_rates[24] += fwd_rates[878];

  //rxn 879
  //sp 25
  sp_rates[24] -= (fwd_rates[879] - rev_rates[760]);
  //sp 4
  sp_rates[3] += (fwd_rates[879] - rev_rates[760]);
  //sp 45
  sp_rates[44] += (fwd_rates[879] - rev_rates[760]);
  //sp 15
  sp_rates[14] -= (fwd_rates[879] - rev_rates[760]);

  //rxn 880
  //sp 26
  sp_rates[25] -= (fwd_rates[880] - rev_rates[761]);
  //sp 11
  sp_rates[10] += (fwd_rates[880] - rev_rates[761]);
  //sp 15
  sp_rates[14] += (fwd_rates[880] - rev_rates[761]);
  //sp 16
  sp_rates[15] -= (fwd_rates[880] - rev_rates[761]);

  //rxn 881
  //sp 3
  sp_rates[2] += (fwd_rates[881] - rev_rates[762]);
  //sp 11
  sp_rates[10] += (fwd_rates[881] - rev_rates[762]);
  //sp 6
  sp_rates[5] -= (fwd_rates[881] - rev_rates[762]);
  //sp 16
  sp_rates[15] -= (fwd_rates[881] - rev_rates[762]);

  //rxn 882
  //sp 25
  sp_rates[24] += (fwd_rates[882] - rev_rates[763]);
  //sp 4
  sp_rates[3] += (fwd_rates[882] - rev_rates[763]);
  //sp 16
  sp_rates[15] -= (fwd_rates[882] - rev_rates[763]);
  //sp 8
  sp_rates[7] -= (fwd_rates[882] - rev_rates[763]);

  //rxn 883
  //sp 25
  sp_rates[24] += (fwd_rates[883] - rev_rates[764]);
  //sp 11
  sp_rates[10] += (fwd_rates[883] - rev_rates[764]);
  //sp 12
  sp_rates[11] -= (fwd_rates[883] - rev_rates[764]);
  //sp 16
  sp_rates[15] -= (fwd_rates[883] - rev_rates[764]);

  //rxn 884
  //sp 17
  sp_rates[16] -= (fwd_rates[884] - rev_rates[765]);
  //sp 11
  sp_rates[10] += (fwd_rates[884] - rev_rates[765]);
  //sp 4
  sp_rates[3] += 2.0 * (fwd_rates[884] - rev_rates[765]);
  //sp 6
  sp_rates[5] -= (fwd_rates[884] - rev_rates[765]);

  //rxn 885
  //sp 81
  sp_rates[80] -= (fwd_rates[885] - rev_rates[766]);
  //sp 26
  sp_rates[25] += (fwd_rates[885] - rev_rates[766]);
  //sp 36
  sp_rates[35] += (fwd_rates[885] - rev_rates[766]);
  //sp 8
  sp_rates[7] -= (fwd_rates[885] - rev_rates[766]);

  //rxn 886
  //sp 19
  sp_rates[18] -= (fwd_rates[886] - rev_rates[767]);
  //sp 4
  sp_rates[3] += (fwd_rates[886] - rev_rates[767]);
  //sp 35
  sp_rates[34] += (fwd_rates[886] - rev_rates[767]);
  //sp 15
  sp_rates[14] -= (fwd_rates[886] - rev_rates[767]);

  //rxn 887
  //sp 18
  sp_rates[17] -= (fwd_rates[887] - rev_rates[768]);
  //sp 4
  sp_rates[3] += (fwd_rates[887] - rev_rates[768]);
  //sp 15
  sp_rates[14] -= (fwd_rates[887] - rev_rates[768]);
  //sp 36
  sp_rates[35] += (fwd_rates[887] - rev_rates[768]);

  //rxn 888
  //sp 17
  sp_rates[16] -= (fwd_rates[888] - rev_rates[769]);
  //sp 14
  sp_rates[13] -= (fwd_rates[888] - rev_rates[769]);
  //sp 15
  sp_rates[14] += 2.0 * (fwd_rates[888] - rev_rates[769]);

  //rxn 889
  //sp 35
  sp_rates[34] += (fwd_rates[889] - rev_rates[770]);
  //sp 36
  sp_rates[35] -= 2.0 * (fwd_rates[889] - rev_rates[770]);
  //sp 37
  sp_rates[36] += (fwd_rates[889] - rev_rates[770]);

  //rxn 890
  //sp 35
  sp_rates[34] += (fwd_rates[890] - rev_rates[771]);
  //sp 36
  sp_rates[35] -= (fwd_rates[890] - rev_rates[771]);
  //sp 30
  sp_rates[29] -= (fwd_rates[890] - rev_rates[771]);
  //sp 31
  sp_rates[30] += (fwd_rates[890] - rev_rates[771]);

  //rxn 891
  //sp 34
  sp_rates[33] += fwd_rates[891];
  //sp 3
  sp_rates[2] += fwd_rates[891];
  //sp 15
  sp_rates[14] -= 2.0 * fwd_rates[891];

  //rxn 892
  //sp 36
  sp_rates[35] += fwd_rates[892];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[892];
  //sp 16
  sp_rates[15] -= 2.0 * fwd_rates[892];

  //rxn 893
  //sp 36
  sp_rates[35] += (fwd_rates[893] - rev_rates[772]);
  //sp 19
  sp_rates[18] -= (fwd_rates[893] - rev_rates[772]);
  //sp 4
  sp_rates[3] += (fwd_rates[893] - rev_rates[772]);
  //sp 16
  sp_rates[15] -= (fwd_rates[893] - rev_rates[772]);

  //rxn 894
  //sp 18
  sp_rates[17] -= (fwd_rates[894] - rev_rates[773]);
  //sp 4
  sp_rates[3] += (fwd_rates[894] - rev_rates[773]);
  //sp 37
  sp_rates[36] += (fwd_rates[894] - rev_rates[773]);
  //sp 16
  sp_rates[15] -= (fwd_rates[894] - rev_rates[773]);

  //rxn 895
  //sp 4
  sp_rates[3] += (fwd_rates[895] - rev_rates[774]);
  //sp 81
  sp_rates[80] += (fwd_rates[895] - rev_rates[774]);
  //sp 19
  sp_rates[18] -= (fwd_rates[895] - rev_rates[774]);
  //sp 36
  sp_rates[35] -= (fwd_rates[895] - rev_rates[774]);

  //rxn 896
  //sp 35
  sp_rates[34] -= fwd_rates[896];
  //sp 102
  sp_rates[101] += fwd_rates[896];
  //sp 62
  sp_rates[61] -= fwd_rates[896];
  //sp 15
  sp_rates[14] += fwd_rates[896];

  //rxn 897
  //sp 34
  sp_rates[33] += fwd_rates[897];
  //sp 59
  sp_rates[58] += fwd_rates[897];
  //sp 60
  sp_rates[59] -= fwd_rates[897];
  //sp 31
  sp_rates[30] -= fwd_rates[897];

  //rxn 898
  //sp 62
  sp_rates[61] += fwd_rates[898];
  //sp 60
  sp_rates[59] -= fwd_rates[898];
  //sp 14
  sp_rates[13] += fwd_rates[898];
  //sp 15
  sp_rates[14] -= fwd_rates[898];

  //rxn 899
  //sp 99
  sp_rates[98] -= fwd_rates[899];
  //sp 102
  sp_rates[101] += fwd_rates[899];
  //sp 30
  sp_rates[29] += fwd_rates[899];
  //sp 31
  sp_rates[30] -= fwd_rates[899];

  //rxn 900
  //sp 63
  sp_rates[62] += (fwd_rates[900] - rev_rates[775]);
  //sp 36
  sp_rates[35] -= (fwd_rates[900] - rev_rates[775]);
  //sp 37
  sp_rates[36] += (fwd_rates[900] - rev_rates[775]);
  //sp 79
  sp_rates[78] -= (fwd_rates[900] - rev_rates[775]);

  //rxn 901
  //sp 35
  sp_rates[34] += (fwd_rates[901] - rev_rates[776]);
  //sp 36
  sp_rates[35] -= (fwd_rates[901] - rev_rates[776]);
  //sp 62
  sp_rates[61] -= (fwd_rates[901] - rev_rates[776]);
  //sp 63
  sp_rates[62] += (fwd_rates[901] - rev_rates[776]);

  //rxn 902
  //sp 35
  sp_rates[34] -= fwd_rates[902];
  //sp 36
  sp_rates[35] += fwd_rates[902];
  //sp 94
  sp_rates[93] += fwd_rates[902];
  //sp 95
  sp_rates[94] -= fwd_rates[902];

  //rxn 903
  //sp 36
  sp_rates[35] += fwd_rates[903];
  //sp 104
  sp_rates[103] += fwd_rates[903];
  //sp 80
  sp_rates[79] -= 2.0 * fwd_rates[903];

  //rxn 904
  //sp 98
  sp_rates[97] += (fwd_rates[904] - rev_rates[777]);
  //sp 99
  sp_rates[98] -= 2.0 * (fwd_rates[904] - rev_rates[777]);
  //sp 102
  sp_rates[101] += (fwd_rates[904] - rev_rates[777]);

  //rxn 905
  //sp 99
  sp_rates[98] -= (fwd_rates[905] - rev_rates[778]);
  //sp 102
  sp_rates[101] += (fwd_rates[905] - rev_rates[778]);
  //sp 94
  sp_rates[93] += (fwd_rates[905] - rev_rates[778]);
  //sp 95
  sp_rates[94] -= (fwd_rates[905] - rev_rates[778]);

  //rxn 906
  //sp 34
  sp_rates[33] += (fwd_rates[906] - rev_rates[779]);
  //sp 35
  sp_rates[34] -= (fwd_rates[906] - rev_rates[779]);
  //sp 102
  sp_rates[101] += (fwd_rates[906] - rev_rates[779]);
  //sp 99
  sp_rates[98] -= (fwd_rates[906] - rev_rates[779]);

  //rxn 907
  //sp 34
  sp_rates[33] += (fwd_rates[907] - rev_rates[780]);
  //sp 37
  sp_rates[36] += (fwd_rates[907] - rev_rates[780]);
  //sp 103
  sp_rates[102] -= (fwd_rates[907] - rev_rates[780]);

  //rxn 908
  //sp 102
  sp_rates[101] += (fwd_rates[908] - rev_rates[781]);
  //sp 15
  sp_rates[14] -= (fwd_rates[908] - rev_rates[781]);
  //sp 80
  sp_rates[79] -= (fwd_rates[908] - rev_rates[781]);

  //rxn 909
  //sp 81
  sp_rates[80] += (fwd_rates[909] - rev_rates[782]);
  //sp 4
  sp_rates[3] += (fwd_rates[909] - rev_rates[782]);
  //sp 80
  sp_rates[79] -= (fwd_rates[909] - rev_rates[782]);

  //rxn 910
  //sp 60
  sp_rates[59] -= (fwd_rates[910] - rev_rates[783]);
  //sp 61
  sp_rates[60] += (fwd_rates[910] - rev_rates[783]);

  //rxn 911
  //sp 63
  sp_rates[62] += fwd_rates[911];
  //sp 34
  sp_rates[33] += fwd_rates[911];
  //sp 102
  sp_rates[101] -= fwd_rates[911];
  //sp 15
  sp_rates[14] -= fwd_rates[911];

  //rxn 912
  //sp 100
  sp_rates[99] += fwd_rates[912];
  //sp 36
  sp_rates[35] -= fwd_rates[912];
  //sp 31
  sp_rates[30] -= fwd_rates[912];

  //rxn 913
  //sp 34
  sp_rates[33] += (fwd_rates[913] - rev_rates[784]);
  //sp 61
  sp_rates[60] += (fwd_rates[913] - rev_rates[784]);
  //sp 62
  sp_rates[61] -= (fwd_rates[913] - rev_rates[784]);
  //sp 31
  sp_rates[30] -= (fwd_rates[913] - rev_rates[784]);

  //rxn 914
  //sp 97
  sp_rates[96] += (fwd_rates[914] - rev_rates[785]);
  //sp 62
  sp_rates[61] -= (fwd_rates[914] - rev_rates[785]);
  //sp 31
  sp_rates[30] -= (fwd_rates[914] - rev_rates[785]);
  //sp 15
  sp_rates[14] += (fwd_rates[914] - rev_rates[785]);

  //rxn 915
  //sp 34
  sp_rates[33] += fwd_rates[915];
  //sp 99
  sp_rates[98] += fwd_rates[915];
  //sp 102
  sp_rates[101] -= fwd_rates[915];
  //sp 31
  sp_rates[30] -= fwd_rates[915];

  //rxn 916
  //sp 34
  sp_rates[33] -= (fwd_rates[916] - rev_rates[786]);
  //sp 60
  sp_rates[59] -= (fwd_rates[916] - rev_rates[786]);
  //sp 62
  sp_rates[61] += (fwd_rates[916] - rev_rates[786]);
  //sp 31
  sp_rates[30] += (fwd_rates[916] - rev_rates[786]);

  //rxn 917
  //sp 60
  sp_rates[59] -= (fwd_rates[917] - rev_rates[787]);
  //sp 61
  sp_rates[60] += (fwd_rates[917] - rev_rates[787]);

  //rxn 918
  //sp 34
  sp_rates[33] -= fwd_rates[918];
  //sp 102
  sp_rates[101] += fwd_rates[918];
  //sp 15
  sp_rates[14] += fwd_rates[918];
  //sp 64
  sp_rates[63] -= fwd_rates[918];

  //rxn 919
  //sp 34
  sp_rates[33] -= fwd_rates[919];
  //sp 99
  sp_rates[98] -= fwd_rates[919];
  //sp 102
  sp_rates[101] += 0.15 * fwd_rates[919];
  //sp 31
  sp_rates[30] += 0.15 * fwd_rates[919];
  //sp 62
  sp_rates[61] += 0.85 * fwd_rates[919];
  //sp 63
  sp_rates[62] += 0.85 * fwd_rates[919];

  //rxn 920
  //sp 34
  sp_rates[33] -= fwd_rates[920];
  //sp 100
  sp_rates[99] -= fwd_rates[920];
  //sp 102
  sp_rates[101] += 0.8 * fwd_rates[920];
  //sp 31
  sp_rates[30] += 0.8 * fwd_rates[920];
  //sp 62
  sp_rates[61] += 0.2 * fwd_rates[920];
  //sp 63
  sp_rates[62] += 0.1 * fwd_rates[920];
  //sp 64
  sp_rates[63] += 0.1 * fwd_rates[920];

  //rxn 921
  //sp 102
  sp_rates[101] += fwd_rates[921];
  //sp 100
  sp_rates[99] -= fwd_rates[921];
  //sp 61
  sp_rates[60] += fwd_rates[921];
  //sp 62
  sp_rates[61] -= fwd_rates[921];

  //rxn 922
  //sp 99
  sp_rates[98] += (fwd_rates[922] - rev_rates[788]);
  //sp 102
  sp_rates[101] -= 2.0 * (fwd_rates[922] - rev_rates[788]);
  //sp 103
  sp_rates[102] += (fwd_rates[922] - rev_rates[788]);

  //rxn 923
  //sp 61
  sp_rates[60] += (fwd_rates[923] - rev_rates[789]);
  //sp 62
  sp_rates[61] -= 2.0 * (fwd_rates[923] - rev_rates[789]);
  //sp 63
  sp_rates[62] += (fwd_rates[923] - rev_rates[789]);

  //rxn 924
  //sp 60
  sp_rates[59] += (fwd_rates[924] - rev_rates[790]);
  //sp 62
  sp_rates[61] -= 2.0 * (fwd_rates[924] - rev_rates[790]);
  //sp 63
  sp_rates[62] += (fwd_rates[924] - rev_rates[790]);

  //rxn 925
  //sp 34
  sp_rates[33] += (fwd_rates[925] - rev_rates[791]);
  //sp 3
  sp_rates[2] += (fwd_rates[925] - rev_rates[791]);
  //sp 30
  sp_rates[29] -= (fwd_rates[925] - rev_rates[791]);

  //rxn 926
  //sp 3
  sp_rates[2] += (fwd_rates[926] - rev_rates[792]);
  //sp 62
  sp_rates[61] -= (fwd_rates[926] - rev_rates[792]);
  //sp 79
  sp_rates[78] += (fwd_rates[926] - rev_rates[792]);

  //rxn 927
  //sp 3
  sp_rates[2] += (fwd_rates[927] - rev_rates[793]);
  //sp 78
  sp_rates[77] += (fwd_rates[927] - rev_rates[793]);
  //sp 62
  sp_rates[61] -= (fwd_rates[927] - rev_rates[793]);

  //rxn 928
  //sp 59
  sp_rates[58] -= (fwd_rates[928] - rev_rates[794]);
  //sp 3
  sp_rates[2] += (fwd_rates[928] - rev_rates[794]);
  //sp 62
  sp_rates[61] += (fwd_rates[928] - rev_rates[794]);

  //rxn 929
  //sp 106
  sp_rates[105] += (fwd_rates[929] - rev_rates[795]);
  //sp 3
  sp_rates[2] += (fwd_rates[929] - rev_rates[795]);
  //sp 104
  sp_rates[103] -= (fwd_rates[929] - rev_rates[795]);

  //rxn 930
  //sp 34
  sp_rates[33] -= (fwd_rates[930] - rev_rates[796]);
  //sp 36
  sp_rates[35] -= (fwd_rates[930] - rev_rates[796]);
  //sp 102
  sp_rates[101] += (fwd_rates[930] - rev_rates[796]);

  //rxn 931
  //sp 97
  sp_rates[96] += (fwd_rates[931] - rev_rates[797]);
  //sp 34
  sp_rates[33] -= 2.0 * (fwd_rates[931] - rev_rates[797]);

  //rxn 932
  //sp 36
  sp_rates[35] -= 2.0 * (fwd_rates[932] - rev_rates[798]);
  //sp 104
  sp_rates[103] += (fwd_rates[932] - rev_rates[798]);

  //rxn 933
  //sp 106
  sp_rates[105] += (fwd_rates[933] - rev_rates[799]);
  //sp 3
  sp_rates[2] += (fwd_rates[933] - rev_rates[799]);
  //sp 36
  sp_rates[35] -= 2.0 * (fwd_rates[933] - rev_rates[799]);

  //rxn 934
  //sp 34
  sp_rates[33] += (fwd_rates[934] - rev_rates[800]);
  //sp 62
  sp_rates[61] -= 2.0 * (fwd_rates[934] - rev_rates[800]);
  //sp 98
  sp_rates[97] += (fwd_rates[934] - rev_rates[800]);

  //rxn 935
  //sp 81
  sp_rates[80] -= 2.0 * (fwd_rates[935] - rev_rates[801]);
  //sp 106
  sp_rates[105] += (fwd_rates[935] - rev_rates[801]);
  //sp 36
  sp_rates[35] += (fwd_rates[935] - rev_rates[801]);

  //rxn 936
  //sp 34
  sp_rates[33] -= fwd_rates[936];
  //sp 107
  sp_rates[106] += fwd_rates[936];
  //sp 106
  sp_rates[105] -= fwd_rates[936];

  //rxn 937
  //sp 4
  sp_rates[3] -= (fwd_rates[937] - rev_rates[802]);
  //sp 79
  sp_rates[78] += (fwd_rates[937] - rev_rates[802]);
  //sp 80
  sp_rates[79] -= (fwd_rates[937] - rev_rates[802]);

  //rxn 938
  //sp 4
  sp_rates[3] -= (fwd_rates[938] - rev_rates[803]);
  //sp 78
  sp_rates[77] += (fwd_rates[938] - rev_rates[803]);
  //sp 80
  sp_rates[79] -= (fwd_rates[938] - rev_rates[803]);

  //rxn 939
  //sp 34
  sp_rates[33] += (fwd_rates[939] - rev_rates[804]);
  //sp 4
  sp_rates[3] -= (fwd_rates[939] - rev_rates[804]);
  //sp 103
  sp_rates[102] -= (fwd_rates[939] - rev_rates[804]);
  //sp 36
  sp_rates[35] += (fwd_rates[939] - rev_rates[804]);

  //rxn 940
  //sp 35
  sp_rates[34] -= (fwd_rates[940] - rev_rates[805]);
  //sp 37
  sp_rates[36] -= (fwd_rates[940] - rev_rates[805]);
  //sp 104
  sp_rates[103] += (fwd_rates[940] - rev_rates[805]);

  //rxn 941
  //sp 36
  sp_rates[35] -= (fwd_rates[941] - rev_rates[806]);
  //sp 37
  sp_rates[36] += (fwd_rates[941] - rev_rates[806]);
  //sp 103
  sp_rates[102] += (fwd_rates[941] - rev_rates[806]);
  //sp 104
  sp_rates[103] -= (fwd_rates[941] - rev_rates[806]);

  //rxn 942
  //sp 105
  sp_rates[104] += (fwd_rates[942] - rev_rates[807]);
  //sp 35
  sp_rates[34] += (fwd_rates[942] - rev_rates[807]);
  //sp 36
  sp_rates[35] -= (fwd_rates[942] - rev_rates[807]);
  //sp 104
  sp_rates[103] -= (fwd_rates[942] - rev_rates[807]);

  //rxn 943
  //sp 105
  sp_rates[104] += (fwd_rates[943] - rev_rates[808]);
  //sp 103
  sp_rates[102] += (fwd_rates[943] - rev_rates[808]);
  //sp 104
  sp_rates[103] -= 2.0 * (fwd_rates[943] - rev_rates[808]);

  //rxn 944
  //sp 17
  sp_rates[16] -= (fwd_rates[944] - rev_rates[809]);
  //sp 4
  sp_rates[3] += (fwd_rates[944] - rev_rates[809]);
  //sp 104
  sp_rates[103] += (fwd_rates[944] - rev_rates[809]);
  //sp 80
  sp_rates[79] -= (fwd_rates[944] - rev_rates[809]);

  //rxn 945
  //sp 81
  sp_rates[80] -= (fwd_rates[945] - rev_rates[810]);
  //sp 4
  sp_rates[3] += (fwd_rates[945] - rev_rates[810]);
  //sp 15
  sp_rates[14] -= (fwd_rates[945] - rev_rates[810]);
  //sp 104
  sp_rates[103] += (fwd_rates[945] - rev_rates[810]);

  //rxn 946
  //sp 99
  sp_rates[98] -= fwd_rates[946];
  //sp 104
  sp_rates[103] += fwd_rates[946];
  //sp 62
  sp_rates[61] += fwd_rates[946];
  //sp 80
  sp_rates[79] -= fwd_rates[946];

  //rxn 947
  //sp 104
  sp_rates[103] += fwd_rates[947];
  //sp 62
  sp_rates[61] += fwd_rates[947];
  //sp 95
  sp_rates[94] -= fwd_rates[947];
  //sp 80
  sp_rates[79] -= fwd_rates[947];

  //rxn 948
  //sp 4
  sp_rates[3] += (fwd_rates[948] - rev_rates[811]);
  //sp 36
  sp_rates[35] -= 2.0 * (fwd_rates[948] - rev_rates[811]);
  //sp 105
  sp_rates[104] += (fwd_rates[948] - rev_rates[811]);

  //rxn 949
  //sp 81
  sp_rates[80] -= (fwd_rates[949] - rev_rates[812]);
  //sp 106
  sp_rates[105] += (fwd_rates[949] - rev_rates[812]);
  //sp 19
  sp_rates[18] -= (fwd_rates[949] - rev_rates[812]);
  //sp 4
  sp_rates[3] += (fwd_rates[949] - rev_rates[812]);

  //rxn 950
  //sp 106
  sp_rates[105] += (fwd_rates[950] - rev_rates[813]);
  //sp 19
  sp_rates[18] -= (fwd_rates[950] - rev_rates[813]);
  //sp 4
  sp_rates[3] += 2.0 * (fwd_rates[950] - rev_rates[813]);
  //sp 80
  sp_rates[79] -= (fwd_rates[950] - rev_rates[813]);

  //rxn 951
  //sp 105
  sp_rates[104] += (fwd_rates[951] - rev_rates[814]);
  //sp 106
  sp_rates[105] -= (fwd_rates[951] - rev_rates[814]);
  //sp 36
  sp_rates[35] -= (fwd_rates[951] - rev_rates[814]);
  //sp 37
  sp_rates[36] += (fwd_rates[951] - rev_rates[814]);

  //rxn 952
  //sp 106
  sp_rates[105] += fwd_rates[952];
  //sp 37
  sp_rates[36] -= fwd_rates[952];
  //sp 78
  sp_rates[77] -= fwd_rates[952];
  //sp 15
  sp_rates[14] += fwd_rates[952];

  //rxn 953
  //sp 34
  sp_rates[33] += fwd_rates[953];
  //sp 79
  sp_rates[78] -= 2.0 * fwd_rates[953];
  //sp 104
  sp_rates[103] += fwd_rates[953];

  //rxn 954
  //sp 34
  sp_rates[33] += fwd_rates[954];
  //sp 78
  sp_rates[77] -= 2.0 * fwd_rates[954];
  //sp 104
  sp_rates[103] += fwd_rates[954];

  //rxn 955
  //sp 107
  sp_rates[106] += fwd_rates[955];
  //sp 36
  sp_rates[35] += fwd_rates[955];
  //sp 104
  sp_rates[103] -= 2.0 * fwd_rates[955];

  //rxn 956
  //sp 105
  sp_rates[104] -= (fwd_rates[956] - rev_rates[815]);
  //sp 37
  sp_rates[36] += (fwd_rates[956] - rev_rates[815]);
  //sp 78
  sp_rates[77] += (fwd_rates[956] - rev_rates[815]);
  //sp 16
  sp_rates[15] -= (fwd_rates[956] - rev_rates[815]);

  //rxn 957
  //sp 81
  sp_rates[80] += (fwd_rates[957] - rev_rates[816]);
  //sp 106
  sp_rates[105] -= (fwd_rates[957] - rev_rates[816]);
  //sp 11
  sp_rates[10] += (fwd_rates[957] - rev_rates[816]);
  //sp 6
  sp_rates[5] -= (fwd_rates[957] - rev_rates[816]);

  //rxn 958
  //sp 81
  sp_rates[80] += fwd_rates[958];
  //sp 106
  sp_rates[105] -= fwd_rates[958];
  //sp 26
  sp_rates[25] += fwd_rates[958];
  //sp 8
  sp_rates[7] -= fwd_rates[958];

  //rxn 959
  //sp 6
  sp_rates[5] -= fwd_rates[959];
  //sp 104
  sp_rates[103] -= fwd_rates[959];
  //sp 11
  sp_rates[10] += fwd_rates[959];
  //sp 78
  sp_rates[77] += 0.5 * fwd_rates[959];
  //sp 79
  sp_rates[78] += 0.5 * fwd_rates[959];

  //rxn 960
  //sp 11
  sp_rates[10] += fwd_rates[960];
  //sp 104
  sp_rates[103] -= fwd_rates[960];
  //sp 63
  sp_rates[62] += fwd_rates[960];
  //sp 8
  sp_rates[7] -= fwd_rates[960];

  //rxn 961
  //sp 35
  sp_rates[34] += 0.4 * fwd_rates[961];
  //sp 5
  sp_rates[4] -= fwd_rates[961];
  //sp 103
  sp_rates[102] -= fwd_rates[961];
  //sp 8
  sp_rates[7] += 0.6 * fwd_rates[961];
  //sp 11
  sp_rates[10] += fwd_rates[961];
  //sp 78
  sp_rates[77] += 0.3 * fwd_rates[961];
  //sp 79
  sp_rates[78] += 0.3 * fwd_rates[961];
  //sp 25
  sp_rates[24] += 0.4 * fwd_rates[961];

  //rxn 962
  //sp 8
  sp_rates[7] -= fwd_rates[962];
  //sp 105
  sp_rates[104] -= fwd_rates[962];
  //sp 11
  sp_rates[10] += fwd_rates[962];
  //sp 78
  sp_rates[77] += 0.5 * fwd_rates[962];
  //sp 79
  sp_rates[78] += 0.5 * fwd_rates[962];

  //rxn 963
  //sp 4
  sp_rates[3] -= (fwd_rates[963] - rev_rates[817]);
  //sp 45
  sp_rates[44] += (fwd_rates[963] - rev_rates[817]);
  //sp 46
  sp_rates[45] -= (fwd_rates[963] - rev_rates[817]);

  //rxn 964
  //sp 26
  sp_rates[25] += fwd_rates[964];
  //sp 46
  sp_rates[45] += fwd_rates[964];
  //sp 5
  sp_rates[4] -= fwd_rates[964];
  //sp 78
  sp_rates[77] -= fwd_rates[964];

  //rxn 965
  //sp 47
  sp_rates[46] += (fwd_rates[965] - rev_rates[818]);
  //sp 10
  sp_rates[9] += (fwd_rates[965] - rev_rates[818]);
  //sp 5
  sp_rates[4] -= (fwd_rates[965] - rev_rates[818]);
  //sp 45
  sp_rates[44] -= (fwd_rates[965] - rev_rates[818]);

  //rxn 966
  //sp 3
  sp_rates[2] += (fwd_rates[966] - rev_rates[819]);
  //sp 4
  sp_rates[3] -= (fwd_rates[966] - rev_rates[819]);
  //sp 45
  sp_rates[44] -= (fwd_rates[966] - rev_rates[819]);
  //sp 46
  sp_rates[45] += (fwd_rates[966] - rev_rates[819]);

  //rxn 967
  //sp 82
  sp_rates[81] -= (fwd_rates[967] - rev_rates[820]) * pres_mod[45];
  //sp 26
  sp_rates[25] += (fwd_rates[967] - rev_rates[820]) * pres_mod[45];
  //sp 31
  sp_rates[30] += (fwd_rates[967] - rev_rates[820]) * pres_mod[45];

  //rxn 968
  //sp 47
  sp_rates[46] += (fwd_rates[968] - rev_rates[821]) * pres_mod[46];
  //sp 82
  sp_rates[81] -= (fwd_rates[968] - rev_rates[821]) * pres_mod[46];
  //sp 15
  sp_rates[14] += (fwd_rates[968] - rev_rates[821]) * pres_mod[46];

  //rxn 969
  //sp 82
  sp_rates[81] -= fwd_rates[969];
  //sp 11
  sp_rates[10] += fwd_rates[969];
  //sp 4
  sp_rates[3] += fwd_rates[969];
  //sp 31
  sp_rates[30] += fwd_rates[969];

  //rxn 970
  //sp 82
  sp_rates[81] -= (fwd_rates[970] - rev_rates[822]);
  //sp 11
  sp_rates[10] += (fwd_rates[970] - rev_rates[822]);
  //sp 30
  sp_rates[29] += (fwd_rates[970] - rev_rates[822]);

  //rxn 971
  //sp 82
  sp_rates[81] -= fwd_rates[971];
  //sp 84
  sp_rates[83] = fwd_rates[971];
  //sp 7
  sp_rates[6] += fwd_rates[971];
  //sp 8
  sp_rates[7] -= fwd_rates[971];

  //rxn 972
  //sp 82
  sp_rates[81] -= fwd_rates[972];
  //sp 10
  sp_rates[9] -= fwd_rates[972];
  //sp 84
  sp_rates[83] += fwd_rates[972];
  //sp 9
  sp_rates[8] += fwd_rates[972];

  //rxn 973
  //sp 34
  sp_rates[33] += 0.25 * fwd_rates[973];
  //sp 4
  sp_rates[3] += 0.75 * fwd_rates[973];
  //sp 84
  sp_rates[83] -= fwd_rates[973];
  //sp 85
  sp_rates[84] += 0.75 * fwd_rates[973];
  //sp 26
  sp_rates[25] += 0.25 * fwd_rates[973];

  //rxn 974
  //sp 10
  sp_rates[9] += fwd_rates[974];
  //sp 84
  sp_rates[83] -= fwd_rates[974];
  //sp 5
  sp_rates[4] -= fwd_rates[974];
  //sp 85
  sp_rates[84] += fwd_rates[974];

  //rxn 975
  //sp 38
  sp_rates[37] += fwd_rates[975];
  //sp 39
  sp_rates[38] -= fwd_rates[975];
  //sp 11
  sp_rates[10] += fwd_rates[975];
  //sp 82
  sp_rates[81] -= fwd_rates[975];
  //sp 31
  sp_rates[30] += fwd_rates[975];

  //rxn 976
  //sp 3
  sp_rates[2] += fwd_rates[976];
  //sp 4
  sp_rates[3] -= fwd_rates[976];
  //sp 11
  sp_rates[10] += fwd_rates[976];
  //sp 82
  sp_rates[81] -= fwd_rates[976];
  //sp 31
  sp_rates[30] += fwd_rates[976];

  //rxn 977
  //sp 7
  sp_rates[6] += fwd_rates[977];
  //sp 8
  sp_rates[7] -= fwd_rates[977];
  //sp 11
  sp_rates[10] += fwd_rates[977];
  //sp 82
  sp_rates[81] -= fwd_rates[977];
  //sp 31
  sp_rates[30] += fwd_rates[977];

  //rxn 978
  //sp 9
  sp_rates[8] += fwd_rates[978];
  //sp 10
  sp_rates[9] -= fwd_rates[978];
  //sp 11
  sp_rates[10] += fwd_rates[978];
  //sp 82
  sp_rates[81] -= fwd_rates[978];
  //sp 31
  sp_rates[30] += fwd_rates[978];

  //rxn 979
  //sp 11
  sp_rates[10] += fwd_rates[979];
  //sp 14
  sp_rates[13] += fwd_rates[979];
  //sp 15
  sp_rates[14] -= fwd_rates[979];
  //sp 82
  sp_rates[81] -= fwd_rates[979];
  //sp 31
  sp_rates[30] += fwd_rates[979];

  //rxn 980
  //sp 4
  sp_rates[3] -= fwd_rates[980];
  //sp 22
  sp_rates[21] -= fwd_rates[980];
  //sp 7
  sp_rates[6] += fwd_rates[980];
  //sp 15
  sp_rates[14] += fwd_rates[980];

  //rxn 981
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[981];
  //sp 62
  sp_rates[61] -= fwd_rates[981];
  //sp 85
  sp_rates[84] += fwd_rates[981];
  //sp 6
  sp_rates[5] -= fwd_rates[981];

  //rxn 982
  //sp 62
  sp_rates[61] += fwd_rates[982];
  //sp 102
  sp_rates[101] -= fwd_rates[982];
  //sp 11
  sp_rates[10] += fwd_rates[982];
  //sp 6
  sp_rates[5] -= fwd_rates[982];

  //rxn 983
  //sp 25
  sp_rates[24] -= fwd_rates[983];
  //sp 12
  sp_rates[11] += fwd_rates[983];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[983];
  //sp 6
  sp_rates[5] -= fwd_rates[983];

  //rxn 984
  //sp 4
  sp_rates[3] += fwd_rates[984];
  //sp 6
  sp_rates[5] -= fwd_rates[984];
  //sp 12
  sp_rates[11] += fwd_rates[984];
  //sp 45
  sp_rates[44] -= fwd_rates[984];
  //sp 15
  sp_rates[14] += fwd_rates[984];

  //rxn 985
  //sp 3
  sp_rates[2] += fwd_rates[985];
  //sp 4
  sp_rates[3] += fwd_rates[985];
  //sp 8
  sp_rates[7] -= fwd_rates[985];
  //sp 12
  sp_rates[11] += fwd_rates[985];
  //sp 25
  sp_rates[24] -= fwd_rates[985];

  //rxn 986
  //sp 3
  sp_rates[2] += fwd_rates[986];
  //sp 8
  sp_rates[7] -= fwd_rates[986];
  //sp 12
  sp_rates[11] += fwd_rates[986];
  //sp 45
  sp_rates[44] -= fwd_rates[986];
  //sp 15
  sp_rates[14] += fwd_rates[986];

  //rxn 987
  //sp 97
  sp_rates[96] -= fwd_rates[987];
  //sp 8
  sp_rates[7] -= fwd_rates[987];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[987];
  //sp 25
  sp_rates[24] += 0.5 * fwd_rates[987];
  //sp 61
  sp_rates[60] += 0.5 * fwd_rates[987];
  //sp 31
  sp_rates[30] += 0.5 * fwd_rates[987];

  //rxn 988
  //sp 98
  sp_rates[97] -= fwd_rates[988];
  //sp 8
  sp_rates[7] -= fwd_rates[988];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[988];
  //sp 15
  sp_rates[14] += 0.5 * fwd_rates[988];
  //sp 82
  sp_rates[81] += 0.5 * fwd_rates[988];
  //sp 31
  sp_rates[30] += 0.5 * fwd_rates[988];

  //rxn 989
  //sp 25
  sp_rates[24] += fwd_rates[989];
  //sp 10
  sp_rates[9] -= fwd_rates[989];
  //sp 36
  sp_rates[35] -= fwd_rates[989];
  //sp 26
  sp_rates[25] += fwd_rates[989];

  //rxn 990
  //sp 26
  sp_rates[25] -= fwd_rates[990];
  //sp 11
  sp_rates[10] += fwd_rates[990];
  //sp 102
  sp_rates[101] -= fwd_rates[990];
  //sp 99
  sp_rates[98] += fwd_rates[990];

  //rxn 991
  //sp 25
  sp_rates[24] -= (fwd_rates[991] - rev_rates[823]);
  //sp 26
  sp_rates[25] -= (fwd_rates[991] - rev_rates[823]);
  //sp 12
  sp_rates[11] += (fwd_rates[991] - rev_rates[823]);
  //sp 15
  sp_rates[14] += (fwd_rates[991] - rev_rates[823]);

  //rxn 992
  //sp 26
  sp_rates[25] -= fwd_rates[992];
  //sp 12
  sp_rates[11] += fwd_rates[992];
  //sp 45
  sp_rates[44] -= fwd_rates[992];
  //sp 31
  sp_rates[30] += fwd_rates[992];

  //rxn 993
  //sp 35
  sp_rates[34] -= fwd_rates[993];
  //sp 36
  sp_rates[35] += fwd_rates[993];
  //sp 23
  sp_rates[22] += 0.5 * fwd_rates[993];
  //sp 24
  sp_rates[23] += 0.5 * fwd_rates[993];
  //sp 25
  sp_rates[24] -= fwd_rates[993];

  //rxn 994
  //sp 11
  sp_rates[10] -= (fwd_rates[994] - rev_rates[824]);
  //sp 12
  sp_rates[11] += (fwd_rates[994] - rev_rates[824]);
  //sp 23
  sp_rates[22] -= (fwd_rates[994] - rev_rates[824]);
  //sp 15
  sp_rates[14] += (fwd_rates[994] - rev_rates[824]);

  //rxn 995
  //sp 25
  sp_rates[24] += fwd_rates[995];
  //sp 11
  sp_rates[10] += fwd_rates[995];
  //sp 36
  sp_rates[35] -= fwd_rates[995];
  //sp 5
  sp_rates[4] -= fwd_rates[995];

  //rxn 996
  //sp 25
  sp_rates[24] += 2.0 * fwd_rates[996];
  //sp 34
  sp_rates[33] -= fwd_rates[996];
  //sp 5
  sp_rates[4] -= fwd_rates[996];

  //rxn 997
  //sp 26
  sp_rates[25] += 2.0 * fwd_rates[997];
  //sp 36
  sp_rates[35] -= fwd_rates[997];
  //sp 5
  sp_rates[4] -= fwd_rates[997];

  //rxn 998
  //sp 34
  sp_rates[33] -= fwd_rates[998];
  //sp 5
  sp_rates[4] -= fwd_rates[998];
  //sp 23
  sp_rates[22] += fwd_rates[998];
  //sp 26
  sp_rates[25] += fwd_rates[998];

  //rxn 999
  //sp 34
  sp_rates[33] += fwd_rates[999];
  //sp 36
  sp_rates[35] += fwd_rates[999];
  //sp 102
  sp_rates[101] -= fwd_rates[999];

  //rxn 1000
  //sp 10
  sp_rates[9] += (fwd_rates[1000] - rev_rates[825]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1000] - rev_rates[825]);
  //sp 78
  sp_rates[77] += (fwd_rates[1000] - rev_rates[825]);
  //sp 64
  sp_rates[63] -= (fwd_rates[1000] - rev_rates[825]);

  //rxn 1001
  //sp 84
  sp_rates[83] += (fwd_rates[1001] - rev_rates[826]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1001] - rev_rates[826]);
  //sp 6
  sp_rates[5] += (fwd_rates[1001] - rev_rates[826]);
  //sp 64
  sp_rates[63] -= (fwd_rates[1001] - rev_rates[826]);

  //rxn 1002
  //sp 10
  sp_rates[9] += fwd_rates[1002];
  //sp 99
  sp_rates[98] -= fwd_rates[1002];
  //sp 5
  sp_rates[4] -= fwd_rates[1002];
  //sp 102
  sp_rates[101] += fwd_rates[1002];

  //rxn 1003
  //sp 10
  sp_rates[9] += fwd_rates[1003];
  //sp 100
  sp_rates[99] -= fwd_rates[1003];
  //sp 5
  sp_rates[4] -= fwd_rates[1003];
  //sp 102
  sp_rates[101] += fwd_rates[1003];

  //rxn 1004
  //sp 5
  sp_rates[4] -= fwd_rates[1004];
  //sp 8
  sp_rates[7] += fwd_rates[1004];
  //sp 11
  sp_rates[10] += fwd_rates[1004];
  //sp 47
  sp_rates[46] -= fwd_rates[1004];
  //sp 25
  sp_rates[24] += fwd_rates[1004];

  //rxn 1005
  //sp 79
  sp_rates[78] += fwd_rates[1005];
  //sp 95
  sp_rates[94] -= fwd_rates[1005];
  //sp 15
  sp_rates[14] += fwd_rates[1005];

  //rxn 1006
  //sp 6
  sp_rates[5] -= fwd_rates[1006] * pres_mod[47];
  //sp 15
  sp_rates[14] -= fwd_rates[1006] * pres_mod[47];
  //sp 23
  sp_rates[22] += fwd_rates[1006] * pres_mod[47];

  //rxn 1007
  //sp 35
  sp_rates[34] -= fwd_rates[1007];
  //sp 6
  sp_rates[5] -= fwd_rates[1007];
  //sp 47
  sp_rates[46] += fwd_rates[1007];

  //rxn 1008
  //sp 25
  sp_rates[24] += fwd_rates[1008];
  //sp 61
  sp_rates[60] -= fwd_rates[1008];
  //sp 6
  sp_rates[5] -= fwd_rates[1008];
  //sp 31
  sp_rates[30] += fwd_rates[1008];

  //rxn 1009
  //sp 25
  sp_rates[24] += fwd_rates[1009];
  //sp 35
  sp_rates[34] += fwd_rates[1009];
  //sp 6
  sp_rates[5] -= fwd_rates[1009];
  //sp 63
  sp_rates[62] -= fwd_rates[1009];

  //rxn 1010
  //sp 99
  sp_rates[98] -= fwd_rates[1010];
  //sp 6
  sp_rates[5] -= fwd_rates[1010];
  //sp 102
  sp_rates[101] += 0.5 * fwd_rates[1010];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1010];
  //sp 25
  sp_rates[24] += 0.5 * fwd_rates[1010];
  //sp 64
  sp_rates[63] += 0.5 * fwd_rates[1010];

  //rxn 1011
  //sp 25
  sp_rates[24] += fwd_rates[1011];
  //sp 65
  sp_rates[64] += fwd_rates[1011];
  //sp 6
  sp_rates[5] -= fwd_rates[1011];
  //sp 95
  sp_rates[94] -= fwd_rates[1011];

  //rxn 1012
  //sp 25
  sp_rates[24] += fwd_rates[1012];
  //sp 26
  sp_rates[25] += fwd_rates[1012];
  //sp 6
  sp_rates[5] -= fwd_rates[1012];
  //sp 47
  sp_rates[46] -= fwd_rates[1012];

  //rxn 1013
  //sp 36
  sp_rates[35] += 0.7 * fwd_rates[1013];
  //sp 37
  sp_rates[36] -= fwd_rates[1013];
  //sp 6
  sp_rates[5] += 0.7 * fwd_rates[1013];
  //sp 8
  sp_rates[7] -= fwd_rates[1013];
  //sp 48
  sp_rates[47] += 0.3 * fwd_rates[1013];

  //rxn 1014
  //sp 35
  sp_rates[34] -= fwd_rates[1014];
  //sp 45
  sp_rates[44] += fwd_rates[1014];
  //sp 8
  sp_rates[7] -= fwd_rates[1014];

  //rxn 1015
  //sp 25
  sp_rates[24] += fwd_rates[1015];
  //sp 62
  sp_rates[61] += fwd_rates[1015];
  //sp 95
  sp_rates[94] -= fwd_rates[1015];
  //sp 8
  sp_rates[7] -= fwd_rates[1015];

  //rxn 1016
  //sp 34
  sp_rates[33] += fwd_rates[1016];
  //sp 99
  sp_rates[98] -= fwd_rates[1016];
  //sp 45
  sp_rates[44] += fwd_rates[1016];
  //sp 8
  sp_rates[7] -= fwd_rates[1016];

  //rxn 1017
  //sp 34
  sp_rates[33] += 0.5 * fwd_rates[1017];
  //sp 99
  sp_rates[98] -= fwd_rates[1017];
  //sp 3
  sp_rates[2] += 0.75 * fwd_rates[1017];
  //sp 8
  sp_rates[7] -= fwd_rates[1017];
  //sp 11
  sp_rates[10] += fwd_rates[1017];
  //sp 62
  sp_rates[61] += 0.5 * fwd_rates[1017];
  //sp 30
  sp_rates[29] += 0.25 * fwd_rates[1017];

  //rxn 1018
  //sp 10
  sp_rates[9] -= fwd_rates[1018];
  //sp 35
  sp_rates[34] -= fwd_rates[1018];
  //sp 47
  sp_rates[46] += fwd_rates[1018];
  //sp 8
  sp_rates[7] += fwd_rates[1018];

  //rxn 1019
  //sp 65
  sp_rates[64] += fwd_rates[1019];
  //sp 8
  sp_rates[7] += fwd_rates[1019];
  //sp 10
  sp_rates[9] -= fwd_rates[1019];
  //sp 25
  sp_rates[24] += fwd_rates[1019];
  //sp 95
  sp_rates[94] -= fwd_rates[1019];

  //rxn 1020
  //sp 9
  sp_rates[8] += (fwd_rates[1020] - rev_rates[827]);
  //sp 10
  sp_rates[9] -= (fwd_rates[1020] - rev_rates[827]);
  //sp 11
  sp_rates[10] += (fwd_rates[1020] - rev_rates[827]);
  //sp 26
  sp_rates[25] -= (fwd_rates[1020] - rev_rates[827]);

  //rxn 1021
  //sp 25
  sp_rates[24] += fwd_rates[1021];
  //sp 7
  sp_rates[6] += fwd_rates[1021];
  //sp 23
  sp_rates[22] -= fwd_rates[1021];
  //sp 8
  sp_rates[7] -= fwd_rates[1021];

  //rxn 1022
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[1022] - rev_rates[828]);
  //sp 26
  sp_rates[25] -= (fwd_rates[1022] - rev_rates[828]);
  //sp 23
  sp_rates[22] -= (fwd_rates[1022] - rev_rates[828]);

  //rxn 1023
  //sp 26
  sp_rates[25] -= fwd_rates[1023];
  //sp 11
  sp_rates[10] += fwd_rates[1023];
  //sp 22
  sp_rates[21] += fwd_rates[1023];
  //sp 23
  sp_rates[22] -= fwd_rates[1023];

  //rxn 1024
  //sp 25
  sp_rates[24] += fwd_rates[1024];
  //sp 14
  sp_rates[13] += fwd_rates[1024];
  //sp 15
  sp_rates[14] -= fwd_rates[1024];
  //sp 24
  sp_rates[23] -= fwd_rates[1024];

  //rxn 1025
  //sp 4
  sp_rates[3] -= fwd_rates[1025];
  //sp 45
  sp_rates[44] += fwd_rates[1025];
  //sp 47
  sp_rates[46] -= fwd_rates[1025];

  //rxn 1026
  //sp 26
  sp_rates[25] += fwd_rates[1026];
  //sp 24
  sp_rates[23] += fwd_rates[1026];
  //sp 47
  sp_rates[46] -= fwd_rates[1026];
  //sp 8
  sp_rates[7] -= fwd_rates[1026];

  //rxn 1027
  //sp 46
  sp_rates[45] += (fwd_rates[1027] - rev_rates[829]);
  //sp 47
  sp_rates[46] -= (fwd_rates[1027] - rev_rates[829]);

  //rxn 1028
  //sp 36
  sp_rates[35] -= fwd_rates[1028];
  //sp 47
  sp_rates[46] += fwd_rates[1028];
  //sp 8
  sp_rates[7] -= fwd_rates[1028];

  //rxn 1029
  //sp 82
  sp_rates[81] += fwd_rates[1029];
  //sp 62
  sp_rates[61] -= fwd_rates[1029];
  //sp 6
  sp_rates[5] -= fwd_rates[1029];

  //rxn 1030
  //sp 97
  sp_rates[96] -= (fwd_rates[1030] - rev_rates[830]);
  //sp 25
  sp_rates[24] += (fwd_rates[1030] - rev_rates[830]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1030] - rev_rates[830]);
  //sp 82
  sp_rates[81] += (fwd_rates[1030] - rev_rates[830]);

  //rxn 1031
  //sp 98
  sp_rates[97] -= (fwd_rates[1031] - rev_rates[831]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1031] - rev_rates[831]);
  //sp 45
  sp_rates[44] += 2.0 * (fwd_rates[1031] - rev_rates[831]);

  //rxn 1032
  //sp 5
  sp_rates[4] -= fwd_rates[1032];
  //sp 102
  sp_rates[101] -= fwd_rates[1032];
  //sp 11
  sp_rates[10] += fwd_rates[1032];
  //sp 26
  sp_rates[25] += fwd_rates[1032];
  //sp 31
  sp_rates[30] += fwd_rates[1032];

  //rxn 1033
  //sp 25
  sp_rates[24] += fwd_rates[1033];
  //sp 82
  sp_rates[81] -= fwd_rates[1033];
  //sp 4
  sp_rates[3] -= fwd_rates[1033];
  //sp 31
  sp_rates[30] += fwd_rates[1033];

  //rxn 1034
  //sp 34
  sp_rates[33] += 0.05 * fwd_rates[1034];
  //sp 35
  sp_rates[34] += 0.4 * fwd_rates[1034];
  //sp 102
  sp_rates[101] -= fwd_rates[1034];
  //sp 8
  sp_rates[7] -= fwd_rates[1034];
  //sp 45
  sp_rates[44] += 0.4 * fwd_rates[1034];
  //sp 47
  sp_rates[46] += 0.05 * fwd_rates[1034];
  //sp 25
  sp_rates[24] += 0.55 * fwd_rates[1034];
  //sp 64
  sp_rates[63] += 0.55 * fwd_rates[1034];

  //rxn 1035
  //sp 35
  sp_rates[34] += fwd_rates[1035];
  //sp 7
  sp_rates[6] += fwd_rates[1035];
  //sp 8
  sp_rates[7] -= fwd_rates[1035];
  //sp 11
  sp_rates[10] += fwd_rates[1035];
  //sp 85
  sp_rates[84] -= fwd_rates[1035];

  //rxn 1036
  //sp 105
  sp_rates[104] += (fwd_rates[1036] - rev_rates[832]);
  //sp 10
  sp_rates[9] += (fwd_rates[1036] - rev_rates[832]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1036] - rev_rates[832]);
  //sp 104
  sp_rates[103] -= (fwd_rates[1036] - rev_rates[832]);

  //rxn 1037
  //sp 4
  sp_rates[3] -= (fwd_rates[1037] - rev_rates[833]);
  //sp 110
  sp_rates[109] += (fwd_rates[1037] - rev_rates[833]);
  //sp 111
  sp_rates[110] = -(fwd_rates[1037] - rev_rates[833]);

  //rxn 1038
  //sp 105
  sp_rates[104] -= fwd_rates[1038];
  //sp 110
  sp_rates[109] += fwd_rates[1038];
  //sp 15
  sp_rates[14] -= fwd_rates[1038];

  //rxn 1039
  //sp 36
  sp_rates[35] -= (fwd_rates[1039] - rev_rates[834]);
  //sp 110
  sp_rates[109] += (fwd_rates[1039] - rev_rates[834]);
  //sp 79
  sp_rates[78] -= (fwd_rates[1039] - rev_rates[834]);

  //rxn 1040
  //sp 110
  sp_rates[109] += (fwd_rates[1040] - rev_rates[835]);
  //sp 36
  sp_rates[35] -= (fwd_rates[1040] - rev_rates[835]);
  //sp 78
  sp_rates[77] -= (fwd_rates[1040] - rev_rates[835]);

  //rxn 1041
  //sp 34
  sp_rates[33] -= (fwd_rates[1041] - rev_rates[836]);
  //sp 4
  sp_rates[3] += (fwd_rates[1041] - rev_rates[836]);
  //sp 110
  sp_rates[109] += (fwd_rates[1041] - rev_rates[836]);
  //sp 80
  sp_rates[79] -= (fwd_rates[1041] - rev_rates[836]);

  //rxn 1042
  //sp 36
  sp_rates[35] += fwd_rates[1042];
  //sp 4
  sp_rates[3] -= fwd_rates[1042];
  //sp 110
  sp_rates[109] -= fwd_rates[1042];
  //sp 63
  sp_rates[62] += fwd_rates[1042];

  //rxn 1043
  //sp 26
  sp_rates[25] += fwd_rates[1043];
  //sp 102
  sp_rates[101] += fwd_rates[1043];
  //sp 110
  sp_rates[109] -= fwd_rates[1043];
  //sp 8
  sp_rates[7] -= fwd_rates[1043];

  //rxn 1044
  //sp 34
  sp_rates[33] += fwd_rates[1044];
  //sp 36
  sp_rates[35] += fwd_rates[1044];
  //sp 8
  sp_rates[7] -= fwd_rates[1044];
  //sp 110
  sp_rates[109] -= fwd_rates[1044];
  //sp 26
  sp_rates[25] += fwd_rates[1044];

  //rxn 1045
  //sp 36
  sp_rates[35] += fwd_rates[1045];
  //sp 8
  sp_rates[7] -= fwd_rates[1045];
  //sp 110
  sp_rates[109] -= fwd_rates[1045];
  //sp 15
  sp_rates[14] += fwd_rates[1045];
  //sp 48
  sp_rates[47] += fwd_rates[1045];

  //rxn 1046
  //sp 11
  sp_rates[10] += 2.0 * fwd_rates[1046];
  //sp 62
  sp_rates[61] += fwd_rates[1046];
  //sp 5
  sp_rates[4] -= fwd_rates[1046];
  //sp 110
  sp_rates[109] -= fwd_rates[1046];

  //rxn 1047
  //sp 5
  sp_rates[4] -= fwd_rates[1047];
  //sp 11
  sp_rates[10] += fwd_rates[1047];
  //sp 110
  sp_rates[109] -= fwd_rates[1047];
  //sp 26
  sp_rates[25] += fwd_rates[1047];
  //sp 63
  sp_rates[62] += fwd_rates[1047];

  //rxn 1048
  //sp 102
  sp_rates[101] += fwd_rates[1048];
  //sp 110
  sp_rates[109] -= fwd_rates[1048];
  //sp 11
  sp_rates[10] += fwd_rates[1048];
  //sp 6
  sp_rates[5] -= fwd_rates[1048];

  //rxn 1049
  //sp 110
  sp_rates[109] -= fwd_rates[1049];
  //sp 6
  sp_rates[5] -= fwd_rates[1049];
  //sp 79
  sp_rates[78] += fwd_rates[1049];
  //sp 48
  sp_rates[47] += fwd_rates[1049];

  //rxn 1050
  //sp 110
  sp_rates[109] -= fwd_rates[1050];
  //sp 78
  sp_rates[77] += fwd_rates[1050];
  //sp 6
  sp_rates[5] -= fwd_rates[1050];
  //sp 48
  sp_rates[47] += fwd_rates[1050];

  //rxn 1051
  //sp 103
  sp_rates[102] += fwd_rates[1051];
  //sp 8
  sp_rates[7] += fwd_rates[1051];
  //sp 10
  sp_rates[9] -= fwd_rates[1051];
  //sp 11
  sp_rates[10] += fwd_rates[1051];
  //sp 111
  sp_rates[110] -= fwd_rates[1051];

  //rxn 1052
  //sp 34
  sp_rates[33] += fwd_rates[1052];
  //sp 37
  sp_rates[36] += fwd_rates[1052];
  //sp 8
  sp_rates[7] += fwd_rates[1052];
  //sp 10
  sp_rates[9] -= fwd_rates[1052];
  //sp 11
  sp_rates[10] += fwd_rates[1052];
  //sp 111
  sp_rates[110] -= fwd_rates[1052];

  //rxn 1053
  //sp 104
  sp_rates[103] += fwd_rates[1053];
  //sp 8
  sp_rates[7] += fwd_rates[1053];
  //sp 10
  sp_rates[9] -= fwd_rates[1053];
  //sp 111
  sp_rates[110] -= fwd_rates[1053];
  //sp 26
  sp_rates[25] += fwd_rates[1053];

  //rxn 1054
  //sp 34
  sp_rates[33] += fwd_rates[1054];
  //sp 36
  sp_rates[35] += fwd_rates[1054];
  //sp 8
  sp_rates[7] -= fwd_rates[1054];
  //sp 11
  sp_rates[10] += fwd_rates[1054];
  //sp 111
  sp_rates[110] -= fwd_rates[1054];

  //rxn 1055
  //sp 4
  sp_rates[3] += fwd_rates[1055];
  //sp 8
  sp_rates[7] -= fwd_rates[1055];
  //sp 111
  sp_rates[110] -= fwd_rates[1055];
  //sp 80
  sp_rates[79] += fwd_rates[1055];
  //sp 48
  sp_rates[47] += fwd_rates[1055];

  //rxn 1056
  //sp 11
  sp_rates[10] += fwd_rates[1056];
  //sp 102
  sp_rates[101] += fwd_rates[1056];
  //sp 111
  sp_rates[110] -= fwd_rates[1056];
  //sp 8
  sp_rates[7] -= fwd_rates[1056];

  //rxn 1057
  //sp 11
  sp_rates[10] += fwd_rates[1057];
  //sp 103
  sp_rates[102] += fwd_rates[1057];
  //sp 6
  sp_rates[5] -= fwd_rates[1057];
  //sp 111
  sp_rates[110] -= fwd_rates[1057];

  //rxn 1058
  //sp 49
  sp_rates[48] += fwd_rates[1058];
  //sp 79
  sp_rates[78] += fwd_rates[1058];
  //sp 6
  sp_rates[5] -= fwd_rates[1058];
  //sp 111
  sp_rates[110] -= fwd_rates[1058];

  //rxn 1059
  //sp 49
  sp_rates[48] += fwd_rates[1059];
  //sp 78
  sp_rates[77] += fwd_rates[1059];
  //sp 6
  sp_rates[5] -= fwd_rates[1059];
  //sp 111
  sp_rates[110] -= fwd_rates[1059];

  //rxn 1060
  //sp 3
  sp_rates[2] += fwd_rates[1060];
  //sp 4
  sp_rates[3] -= fwd_rates[1060];
  //sp 110
  sp_rates[109] -= fwd_rates[1060];
  //sp 111
  sp_rates[110] += fwd_rates[1060];

  //rxn 1061
  //sp 3
  sp_rates[2] -= fwd_rates[1061];
  //sp 4
  sp_rates[3] += fwd_rates[1061];
  //sp 110
  sp_rates[109] += fwd_rates[1061];
  //sp 111
  sp_rates[110] -= fwd_rates[1061];

  //rxn 1062
  //sp 14
  sp_rates[13] += fwd_rates[1062];
  //sp 111
  sp_rates[110] += fwd_rates[1062];
  //sp 110
  sp_rates[109] -= fwd_rates[1062];
  //sp 15
  sp_rates[14] -= fwd_rates[1062];

  //rxn 1063
  //sp 110
  sp_rates[109] += fwd_rates[1063];
  //sp 14
  sp_rates[13] -= fwd_rates[1063];
  //sp 111
  sp_rates[110] -= fwd_rates[1063];
  //sp 15
  sp_rates[14] += fwd_rates[1063];

  //rxn 1064
  //sp 111
  sp_rates[110] += fwd_rates[1064];
  //sp 110
  sp_rates[109] -= fwd_rates[1064];
  //sp 7
  sp_rates[6] += fwd_rates[1064];
  //sp 8
  sp_rates[7] -= fwd_rates[1064];

  //rxn 1065
  //sp 36
  sp_rates[35] += fwd_rates[1065];
  //sp 111
  sp_rates[110] -= fwd_rates[1065];
  //sp 80
  sp_rates[79] += fwd_rates[1065];

  //rxn 1066
  //sp 106
  sp_rates[105] -= fwd_rates[1066];
  //sp 111
  sp_rates[110] += fwd_rates[1066];
  //sp 15
  sp_rates[14] -= fwd_rates[1066];

  //rxn 1067
  //sp 103
  sp_rates[102] += fwd_rates[1067];
  //sp 110
  sp_rates[109] -= fwd_rates[1067];
  //sp 111
  sp_rates[110] -= fwd_rates[1067];
  //sp 107
  sp_rates[106] += fwd_rates[1067];

  //rxn 1068
  //sp 111
  sp_rates[110] -= fwd_rates[1068];
  //sp 103
  sp_rates[102] += fwd_rates[1068];
  //sp 35
  sp_rates[34] += fwd_rates[1068];
  //sp 15
  sp_rates[14] -= fwd_rates[1068];

  //rxn 1069
  //sp 111
  sp_rates[110] -= fwd_rates[1069];
  //sp 36
  sp_rates[35] += fwd_rates[1069];
  //sp 102
  sp_rates[101] += fwd_rates[1069];
  //sp 15
  sp_rates[14] -= fwd_rates[1069];

  //rxn 1070
  //sp 34
  sp_rates[33] += fwd_rates[1070];
  //sp 111
  sp_rates[110] -= fwd_rates[1070];
  //sp 15
  sp_rates[14] -= fwd_rates[1070];
  //sp 104
  sp_rates[103] += fwd_rates[1070];

  //rxn 1071
  //sp 35
  sp_rates[34] -= fwd_rates[1071];
  //sp 4
  sp_rates[3] += fwd_rates[1071];
  //sp 110
  sp_rates[109] += fwd_rates[1071];
  //sp 79
  sp_rates[78] -= fwd_rates[1071];

  //rxn 1072
  //sp 110
  sp_rates[109] += fwd_rates[1072];
  //sp 35
  sp_rates[34] -= fwd_rates[1072];
  //sp 4
  sp_rates[3] += fwd_rates[1072];
  //sp 78
  sp_rates[77] -= fwd_rates[1072];

  //rxn 1073
  //sp 34
  sp_rates[33] -= fwd_rates[1073];
  //sp 99
  sp_rates[98] -= fwd_rates[1073];
  //sp 3
  sp_rates[2] += fwd_rates[1073];
  //sp 110
  sp_rates[109] += fwd_rates[1073];
  //sp 15
  sp_rates[14] += fwd_rates[1073];

  //rxn 1074
  //sp 34
  sp_rates[33] -= fwd_rates[1074];
  //sp 4
  sp_rates[3] += fwd_rates[1074];
  //sp 110
  sp_rates[109] += fwd_rates[1074];
  //sp 14
  sp_rates[13] += fwd_rates[1074];
  //sp 95
  sp_rates[94] -= fwd_rates[1074];

  //rxn 1075
  //sp 34
  sp_rates[33] -= fwd_rates[1075];
  //sp 79
  sp_rates[78] += 0.33 * fwd_rates[1075];
  //sp 95
  sp_rates[94] -= fwd_rates[1075];
  //sp 61
  sp_rates[60] += 0.33 * fwd_rates[1075];
  //sp 62
  sp_rates[61] += 0.67 * fwd_rates[1075];
  //sp 63
  sp_rates[62] += 0.67 * fwd_rates[1075];

  //rxn 1076
  //sp 79
  sp_rates[78] += fwd_rates[1076];
  //sp 111
  sp_rates[110] -= fwd_rates[1076];
  //sp 103
  sp_rates[102] -= fwd_rates[1076];
  //sp 107
  sp_rates[106] += fwd_rates[1076];

  //rxn 1077
  //sp 63
  sp_rates[62] += fwd_rates[1077];
  //sp 107
  sp_rates[106] += fwd_rates[1077];
  //sp 110
  sp_rates[109] -= fwd_rates[1077];
  //sp 103
  sp_rates[102] -= fwd_rates[1077];

  //rxn 1078
  //sp 34
  sp_rates[33] += (fwd_rates[1078] - rev_rates[837]);
  //sp 35
  sp_rates[34] -= (fwd_rates[1078] - rev_rates[837]);
  //sp 110
  sp_rates[109] -= (fwd_rates[1078] - rev_rates[837]);
  //sp 111
  sp_rates[110] += (fwd_rates[1078] - rev_rates[837]);

  //rxn 1079
  //sp 107
  sp_rates[106] += fwd_rates[1079];
  //sp 37
  sp_rates[36] -= fwd_rates[1079];
  //sp 103
  sp_rates[102] -= fwd_rates[1079];

  //rxn 1080
  //sp 105
  sp_rates[104] -= fwd_rates[1080];
  //sp 35
  sp_rates[34] -= fwd_rates[1080];
  //sp 107
  sp_rates[106] += fwd_rates[1080];

  //rxn 1081
  //sp 105
  sp_rates[104] -= (fwd_rates[1081] - rev_rates[838]);
  //sp 35
  sp_rates[34] -= (fwd_rates[1081] - rev_rates[838]);
  //sp 4
  sp_rates[3] += (fwd_rates[1081] - rev_rates[838]);
  //sp 109
  sp_rates[108] += (fwd_rates[1081] - rev_rates[838]);

  //rxn 1082
  //sp 34
  sp_rates[33] += (fwd_rates[1082] - rev_rates[839]);
  //sp 107
  sp_rates[106] += (fwd_rates[1082] - rev_rates[839]);
  //sp 103
  sp_rates[102] -= 2.0 * (fwd_rates[1082] - rev_rates[839]);

  //rxn 1083
  //sp 34
  sp_rates[33] -= fwd_rates[1083];
  //sp 4
  sp_rates[3] += fwd_rates[1083];
  //sp 107
  sp_rates[106] += fwd_rates[1083];
  //sp 110
  sp_rates[109] -= fwd_rates[1083];
  //sp 15
  sp_rates[14] += fwd_rates[1083];

  //rxn 1084
  //sp 110
  sp_rates[109] -= fwd_rates[1084];
  //sp 62
  sp_rates[61] -= fwd_rates[1084];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[1084];
  //sp 107
  sp_rates[106] += fwd_rates[1084];

  //rxn 1085
  //sp 97
  sp_rates[96] -= fwd_rates[1085];
  //sp 107
  sp_rates[106] += fwd_rates[1085];
  //sp 110
  sp_rates[109] -= fwd_rates[1085];
  //sp 15
  sp_rates[14] += fwd_rates[1085];
  //sp 31
  sp_rates[30] += fwd_rates[1085];

  //rxn 1086
  //sp 98
  sp_rates[97] -= fwd_rates[1086];
  //sp 107
  sp_rates[106] += fwd_rates[1086];
  //sp 110
  sp_rates[109] -= fwd_rates[1086];
  //sp 15
  sp_rates[14] += fwd_rates[1086];
  //sp 31
  sp_rates[30] += fwd_rates[1086];

  //rxn 1087
  //sp 35
  sp_rates[34] += fwd_rates[1087];
  //sp 102
  sp_rates[101] -= fwd_rates[1087];
  //sp 107
  sp_rates[106] += fwd_rates[1087];
  //sp 110
  sp_rates[109] -= fwd_rates[1087];
  //sp 15
  sp_rates[14] += fwd_rates[1087];

  //rxn 1088
  //sp 35
  sp_rates[34] -= fwd_rates[1088];
  //sp 4
  sp_rates[3] += fwd_rates[1088];
  //sp 107
  sp_rates[106] += fwd_rates[1088];
  //sp 104
  sp_rates[103] -= fwd_rates[1088];

  //rxn 1089
  //sp 35
  sp_rates[34] -= fwd_rates[1089];
  //sp 4
  sp_rates[3] += fwd_rates[1089];
  //sp 3
  sp_rates[2] += fwd_rates[1089];
  //sp 102
  sp_rates[101] -= fwd_rates[1089];
  //sp 107
  sp_rates[106] += fwd_rates[1089];

  //rxn 1090
  //sp 35
  sp_rates[34] -= fwd_rates[1090];
  //sp 110
  sp_rates[109] += fwd_rates[1090];
  //sp 102
  sp_rates[101] -= fwd_rates[1090];
  //sp 15
  sp_rates[14] += fwd_rates[1090];

  //rxn 1091
  //sp 107
  sp_rates[106] += fwd_rates[1091];
  //sp 4
  sp_rates[3] += fwd_rates[1091];
  //sp 37
  sp_rates[36] -= fwd_rates[1091];
  //sp 102
  sp_rates[101] -= fwd_rates[1091];

  //rxn 1092
  //sp 107
  sp_rates[106] += (fwd_rates[1092] - rev_rates[840]);
  //sp 4
  sp_rates[3] += (fwd_rates[1092] - rev_rates[840]);
  //sp 78
  sp_rates[77] -= (fwd_rates[1092] - rev_rates[840]);
  //sp 80
  sp_rates[79] -= (fwd_rates[1092] - rev_rates[840]);

  //rxn 1093
  //sp 105
  sp_rates[104] -= fwd_rates[1093];
  //sp 34
  sp_rates[33] -= fwd_rates[1093];
  //sp 107
  sp_rates[106] += fwd_rates[1093];
  //sp 4
  sp_rates[3] += fwd_rates[1093];

  //rxn 1094
  //sp 103
  sp_rates[102] -= fwd_rates[1094];
  //sp 107
  sp_rates[106] += fwd_rates[1094];
  //sp 79
  sp_rates[78] -= fwd_rates[1094];
  //sp 15
  sp_rates[14] += fwd_rates[1094];

  //rxn 1095
  //sp 107
  sp_rates[106] += fwd_rates[1095];
  //sp 78
  sp_rates[77] -= fwd_rates[1095];
  //sp 103
  sp_rates[102] -= fwd_rates[1095];
  //sp 15
  sp_rates[14] += fwd_rates[1095];

  //rxn 1096
  //sp 106
  sp_rates[105] -= fwd_rates[1096];
  //sp 107
  sp_rates[106] += fwd_rates[1096];
  //sp 37
  sp_rates[36] += fwd_rates[1096];
  //sp 103
  sp_rates[102] -= fwd_rates[1096];

  //rxn 1097
  //sp 35
  sp_rates[34] += fwd_rates[1097];
  //sp 107
  sp_rates[106] += fwd_rates[1097];
  //sp 103
  sp_rates[102] -= fwd_rates[1097];
  //sp 104
  sp_rates[103] -= fwd_rates[1097];

  //rxn 1098
  //sp 11
  sp_rates[10] += fwd_rates[1098];
  //sp 109
  sp_rates[108] -= fwd_rates[1098];
  //sp 6
  sp_rates[5] -= fwd_rates[1098];
  //sp 111
  sp_rates[110] += fwd_rates[1098];

  //rxn 1099
  //sp 107
  sp_rates[106] -= (fwd_rates[1099] - rev_rates[841]);
  //sp 109
  sp_rates[108] += (fwd_rates[1099] - rev_rates[841]);
  //sp 86
  sp_rates[85] += (fwd_rates[1099] - rev_rates[841]);
  //sp 87
  sp_rates[86] -= (fwd_rates[1099] - rev_rates[841]);

  //rxn 1100
  //sp 113
  sp_rates[112] = fwd_rates[1100];
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[1100];
  //sp 110
  sp_rates[109] -= 2.0 * fwd_rates[1100];

  //rxn 1101
  //sp 3
  sp_rates[2] += fwd_rates[1101];
  //sp 4
  sp_rates[3] += fwd_rates[1101];
  //sp 110
  sp_rates[109] -= fwd_rates[1101];
  //sp 111
  sp_rates[110] -= fwd_rates[1101];
  //sp 113
  sp_rates[112] += fwd_rates[1101];

  //rxn 1102
  //sp 113
  sp_rates[112] += fwd_rates[1102];
  //sp 107
  sp_rates[106] -= fwd_rates[1102];
  //sp 111
  sp_rates[110] -= fwd_rates[1102];
  //sp 15
  sp_rates[14] += fwd_rates[1102];

  //rxn 1103
  //sp 113
  sp_rates[112] += fwd_rates[1103];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[1103];
  //sp 111
  sp_rates[110] -= 2.0 * fwd_rates[1103];

  //rxn 1104
  //sp 105
  sp_rates[104] -= fwd_rates[1104];
  //sp 107
  sp_rates[106] -= fwd_rates[1104];
  //sp 4
  sp_rates[3] += fwd_rates[1104];
  //sp 113
  sp_rates[112] += fwd_rates[1104];

  //rxn 1105
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[1105];
  //sp 4
  sp_rates[3] += fwd_rates[1105];
  //sp 100
  sp_rates[99] -= fwd_rates[1105];
  //sp 107
  sp_rates[106] -= fwd_rates[1105];
  //sp 113
  sp_rates[112] += fwd_rates[1105];

  //rxn 1106
  //sp 3
  sp_rates[2] += fwd_rates[1106];
  //sp 4
  sp_rates[3] += 2.0 * fwd_rates[1106];
  //sp 102
  sp_rates[101] -= fwd_rates[1106];
  //sp 107
  sp_rates[106] -= fwd_rates[1106];
  //sp 113
  sp_rates[112] += fwd_rates[1106];

  //rxn 1107
  //sp 4
  sp_rates[3] += fwd_rates[1107];
  //sp 107
  sp_rates[106] -= fwd_rates[1107];
  //sp 110
  sp_rates[109] -= fwd_rates[1107];
  //sp 15
  sp_rates[14] += fwd_rates[1107];
  //sp 113
  sp_rates[112] += fwd_rates[1107];

  //rxn 1108
  //sp 113
  sp_rates[112] += fwd_rates[1108];
  //sp 4
  sp_rates[3] += fwd_rates[1108];
  //sp 109
  sp_rates[108] -= fwd_rates[1108];
  //sp 104
  sp_rates[103] -= fwd_rates[1108];

  //rxn 1109
  //sp 113
  sp_rates[112] -= fwd_rates[1109];
  //sp 4
  sp_rates[3] -= fwd_rates[1109];
  //sp 109
  sp_rates[108] += fwd_rates[1109];
  //sp 104
  sp_rates[103] += fwd_rates[1109];

  //rxn 1110
  //sp 111
  sp_rates[110] -= (fwd_rates[1110] - rev_rates[842]);
  //sp 15
  sp_rates[14] -= (fwd_rates[1110] - rev_rates[842]);
  //sp 112
  sp_rates[111] = (fwd_rates[1110] - rev_rates[842]);

  //rxn 1111
  //sp 108
  sp_rates[107] += fwd_rates[1111];
  //sp 111
  sp_rates[110] -= fwd_rates[1111];
  //sp 3
  sp_rates[2] += fwd_rates[1111];
  //sp 15
  sp_rates[14] -= fwd_rates[1111];

  //rxn 1112
  //sp 3
  sp_rates[2] += fwd_rates[1112];
  //sp 108
  sp_rates[107] += fwd_rates[1112];
  //sp 112
  sp_rates[111] -= fwd_rates[1112];

  //rxn 1113
  //sp 4
  sp_rates[3] -= fwd_rates[1113];
  //sp 110
  sp_rates[109] += fwd_rates[1113];
  //sp 15
  sp_rates[14] += fwd_rates[1113];
  //sp 112
  sp_rates[111] -= fwd_rates[1113];

  //rxn 1114
  //sp 78
  sp_rates[77] -= (fwd_rates[1114] - rev_rates[843]);
  //sp 79
  sp_rates[78] -= (fwd_rates[1114] - rev_rates[843]);
  //sp 112
  sp_rates[111] += (fwd_rates[1114] - rev_rates[843]);

  //rxn 1115
  //sp 4
  sp_rates[3] += fwd_rates[1115];
  //sp 36
  sp_rates[35] -= fwd_rates[1115];
  //sp 95
  sp_rates[94] -= fwd_rates[1115];
  //sp 112
  sp_rates[111] += fwd_rates[1115];

  //rxn 1116
  //sp 4
  sp_rates[3] += fwd_rates[1116];
  //sp 99
  sp_rates[98] -= fwd_rates[1116];
  //sp 36
  sp_rates[35] -= fwd_rates[1116];
  //sp 112
  sp_rates[111] += fwd_rates[1116];

  //rxn 1117
  //sp 4
  sp_rates[3] -= (fwd_rates[1117] - rev_rates[844]) * pres_mod[48];
  //sp 62
  sp_rates[61] += (fwd_rates[1117] - rev_rates[844]) * pres_mod[48];
  //sp 63
  sp_rates[62] -= (fwd_rates[1117] - rev_rates[844]) * pres_mod[48];

  //rxn 1118
  //sp 34
  sp_rates[33] += (fwd_rates[1118] - rev_rates[845]);
  //sp 100
  sp_rates[99] += (fwd_rates[1118] - rev_rates[845]);
  //sp 62
  sp_rates[61] -= (fwd_rates[1118] - rev_rates[845]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1118] - rev_rates[845]);

  //rxn 1119
  //sp 3
  sp_rates[2] += fwd_rates[1119];
  //sp 110
  sp_rates[109] += fwd_rates[1119];
  //sp 15
  sp_rates[14] += fwd_rates[1119];
  //sp 62
  sp_rates[61] -= fwd_rates[1119];
  //sp 63
  sp_rates[62] -= fwd_rates[1119];

  //rxn 1120
  //sp 3
  sp_rates[2] += fwd_rates[1120];
  //sp 4
  sp_rates[3] += fwd_rates[1120];
  //sp 112
  sp_rates[111] += fwd_rates[1120];
  //sp 62
  sp_rates[61] -= fwd_rates[1120];
  //sp 63
  sp_rates[62] -= fwd_rates[1120];

  //rxn 1121
  //sp 97
  sp_rates[96] += (fwd_rates[1121] - rev_rates[846]);
  //sp 35
  sp_rates[34] += (fwd_rates[1121] - rev_rates[846]);
  //sp 62
  sp_rates[61] -= (fwd_rates[1121] - rev_rates[846]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1121] - rev_rates[846]);

  //rxn 1122
  //sp 97
  sp_rates[96] -= fwd_rates[1122];
  //sp 3
  sp_rates[2] += fwd_rates[1122];
  //sp 15
  sp_rates[14] += fwd_rates[1122];
  //sp 112
  sp_rates[111] += fwd_rates[1122];
  //sp 63
  sp_rates[62] -= fwd_rates[1122];

  //rxn 1123
  //sp 98
  sp_rates[97] -= fwd_rates[1123];
  //sp 3
  sp_rates[2] += fwd_rates[1123];
  //sp 15
  sp_rates[14] += fwd_rates[1123];
  //sp 112
  sp_rates[111] += fwd_rates[1123];
  //sp 63
  sp_rates[62] -= fwd_rates[1123];

  //rxn 1124
  //sp 110
  sp_rates[109] += (fwd_rates[1124] - rev_rates[847]);
  //sp 102
  sp_rates[101] -= (fwd_rates[1124] - rev_rates[847]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1124] - rev_rates[847]);
  //sp 31
  sp_rates[30] += (fwd_rates[1124] - rev_rates[847]);

  //rxn 1125
  //sp 89
  sp_rates[88] -= (fwd_rates[1125] - rev_rates[848]);
  //sp 34
  sp_rates[33] -= (fwd_rates[1125] - rev_rates[848]);
  //sp 61
  sp_rates[60] += (fwd_rates[1125] - rev_rates[848]);
  //sp 62
  sp_rates[61] += (fwd_rates[1125] - rev_rates[848]);

  //rxn 1126
  //sp 36
  sp_rates[35] += (fwd_rates[1126] - rev_rates[849]);
  //sp 102
  sp_rates[101] += (fwd_rates[1126] - rev_rates[849]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1126] - rev_rates[849]);
  //sp 80
  sp_rates[79] -= (fwd_rates[1126] - rev_rates[849]);

  //rxn 1127
  //sp 99
  sp_rates[98] -= (fwd_rates[1127] - rev_rates[850]);
  //sp 102
  sp_rates[101] += (fwd_rates[1127] - rev_rates[850]);
  //sp 62
  sp_rates[61] += (fwd_rates[1127] - rev_rates[850]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1127] - rev_rates[850]);

  //rxn 1128
  //sp 95
  sp_rates[94] -= (fwd_rates[1128] - rev_rates[851]);
  //sp 94
  sp_rates[93] += (fwd_rates[1128] - rev_rates[851]);
  //sp 63
  sp_rates[62] -= (fwd_rates[1128] - rev_rates[851]);
  //sp 79
  sp_rates[78] += (fwd_rates[1128] - rev_rates[851]);

  //rxn 1129
  //sp 78
  sp_rates[77] -= (fwd_rates[1129] - rev_rates[852]);
  //sp 79
  sp_rates[78] += (fwd_rates[1129] - rev_rates[852]);

  //rxn 1130
  //sp 106
  sp_rates[105] -= fwd_rates[1130];
  //sp 35
  sp_rates[34] -= fwd_rates[1130];
  //sp 109
  sp_rates[108] += fwd_rates[1130];

  //rxn 1131
  //sp 97
  sp_rates[96] -= (fwd_rates[1131] - rev_rates[853]);
  //sp 89
  sp_rates[88] += (fwd_rates[1131] - rev_rates[853]);
  //sp 100
  sp_rates[99] -= (fwd_rates[1131] - rev_rates[853]);
  //sp 102
  sp_rates[101] += (fwd_rates[1131] - rev_rates[853]);

  //rxn 1132
  //sp 90
  sp_rates[89] += (fwd_rates[1132] - rev_rates[854]);
  //sp 98
  sp_rates[97] -= (fwd_rates[1132] - rev_rates[854]);
  //sp 100
  sp_rates[99] -= (fwd_rates[1132] - rev_rates[854]);
  //sp 102
  sp_rates[101] += (fwd_rates[1132] - rev_rates[854]);

  //rxn 1133
  //sp 97
  sp_rates[96] += 0.5 * fwd_rates[1133];
  //sp 98
  sp_rates[97] += 0.5 * fwd_rates[1133];
  //sp 3
  sp_rates[2] += fwd_rates[1133];
  //sp 88
  sp_rates[87] -= fwd_rates[1133];

  //rxn 1134
  //sp 91
  sp_rates[90] -= (fwd_rates[1134] - rev_rates[855]);
  //sp 3
  sp_rates[2] += (fwd_rates[1134] - rev_rates[855]);
  //sp 94
  sp_rates[93] += (fwd_rates[1134] - rev_rates[855]);

  //rxn 1135
  //sp 90
  sp_rates[89] += fwd_rates[1135];
  //sp 26
  sp_rates[25] -= fwd_rates[1135];
  //sp 11
  sp_rates[10] += fwd_rates[1135];
  //sp 98
  sp_rates[97] -= fwd_rates[1135];

  //rxn 1136
  //sp 97
  sp_rates[96] -= fwd_rates[1136];
  //sp 26
  sp_rates[25] -= fwd_rates[1136];
  //sp 11
  sp_rates[10] += fwd_rates[1136];
  //sp 89
  sp_rates[88] += fwd_rates[1136];

  //rxn 1137
  //sp 11
  sp_rates[10] += fwd_rates[1137];
  //sp 26
  sp_rates[25] -= fwd_rates[1137];
  //sp 92
  sp_rates[91] += 0.7 * fwd_rates[1137];
  //sp 93
  sp_rates[92] += 0.3 * fwd_rates[1137];
  //sp 94
  sp_rates[93] -= fwd_rates[1137];

  //rxn 1138
  //sp 25
  sp_rates[24] += fwd_rates[1138];
  //sp 86
  sp_rates[85] += fwd_rates[1138];
  //sp 5
  sp_rates[4] -= fwd_rates[1138];
  //sp 94
  sp_rates[93] -= fwd_rates[1138];

  //rxn 1139
  //sp 10
  sp_rates[9] += (fwd_rates[1139] - rev_rates[856]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1139] - rev_rates[856]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1139] - rev_rates[856]);
  //sp 94
  sp_rates[93] += (fwd_rates[1139] - rev_rates[856]);

  //rxn 1140
  //sp 10
  sp_rates[9] += (fwd_rates[1140] - rev_rates[857]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1140] - rev_rates[857]);
  //sp 94
  sp_rates[93] += (fwd_rates[1140] - rev_rates[857]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1140] - rev_rates[857]);

  //rxn 1141
  //sp 89
  sp_rates[88] -= fwd_rates[1141];
  //sp 25
  sp_rates[24] += fwd_rates[1141];
  //sp 61
  sp_rates[60] += fwd_rates[1141];
  //sp 6
  sp_rates[5] -= fwd_rates[1141];

  //rxn 1142
  //sp 82
  sp_rates[81] += fwd_rates[1142];
  //sp 90
  sp_rates[89] -= fwd_rates[1142];
  //sp 6
  sp_rates[5] -= fwd_rates[1142];
  //sp 15
  sp_rates[14] += fwd_rates[1142];

  //rxn 1143
  //sp 25
  sp_rates[24] += fwd_rates[1143];
  //sp 60
  sp_rates[59] += fwd_rates[1143];
  //sp 92
  sp_rates[91] -= fwd_rates[1143];
  //sp 6
  sp_rates[5] -= fwd_rates[1143];

  //rxn 1144
  //sp 86
  sp_rates[85] += fwd_rates[1144];
  //sp 93
  sp_rates[92] -= fwd_rates[1144];
  //sp 6
  sp_rates[5] -= fwd_rates[1144];
  //sp 15
  sp_rates[14] += fwd_rates[1144];

  //rxn 1145
  //sp 3
  sp_rates[2] -= fwd_rates[1145];
  //sp 4
  sp_rates[3] += fwd_rates[1145];
  //sp 102
  sp_rates[101] += fwd_rates[1145];
  //sp 103
  sp_rates[102] -= fwd_rates[1145];

  //rxn 1146
  //sp 4
  sp_rates[3] += fwd_rates[1146];
  //sp 3
  sp_rates[2] -= fwd_rates[1146];
  //sp 100
  sp_rates[99] -= fwd_rates[1146];
  //sp 97
  sp_rates[96] += fwd_rates[1146];

  //rxn 1147
  //sp 3
  sp_rates[2] -= fwd_rates[1147];
  //sp 4
  sp_rates[3] += fwd_rates[1147];
  //sp 94
  sp_rates[93] += fwd_rates[1147];
  //sp 95
  sp_rates[94] -= fwd_rates[1147];

  //rxn 1148
  //sp 49
  sp_rates[48] -= fwd_rates[1148];
  //sp 3
  sp_rates[2] -= fwd_rates[1148];
  //sp 4
  sp_rates[3] += fwd_rates[1148];
  //sp 48
  sp_rates[47] += fwd_rates[1148];

  //rxn 1149
  //sp 102
  sp_rates[101] += fwd_rates[1149];
  //sp 14
  sp_rates[13] -= fwd_rates[1149];
  //sp 103
  sp_rates[102] -= fwd_rates[1149];
  //sp 15
  sp_rates[14] += fwd_rates[1149];

  //rxn 1150
  //sp 97
  sp_rates[96] += fwd_rates[1150];
  //sp 100
  sp_rates[99] -= fwd_rates[1150];
  //sp 14
  sp_rates[13] -= fwd_rates[1150];
  //sp 15
  sp_rates[14] += fwd_rates[1150];

  //rxn 1151
  //sp 94
  sp_rates[93] += fwd_rates[1151];
  //sp 14
  sp_rates[13] -= fwd_rates[1151];
  //sp 95
  sp_rates[94] -= fwd_rates[1151];
  //sp 15
  sp_rates[14] += fwd_rates[1151];

  //rxn 1152
  //sp 107
  sp_rates[106] += fwd_rates[1152];
  //sp 109
  sp_rates[108] -= fwd_rates[1152];
  //sp 14
  sp_rates[13] -= fwd_rates[1152];
  //sp 15
  sp_rates[14] += fwd_rates[1152];

  //rxn 1153
  //sp 46
  sp_rates[45] -= fwd_rates[1153];
  //sp 45
  sp_rates[44] += fwd_rates[1153];
  //sp 14
  sp_rates[13] -= fwd_rates[1153];
  //sp 15
  sp_rates[14] += fwd_rates[1153];

  //rxn 1154
  //sp 45
  sp_rates[44] += fwd_rates[1154];
  //sp 14
  sp_rates[13] -= fwd_rates[1154];
  //sp 47
  sp_rates[46] -= fwd_rates[1154];
  //sp 15
  sp_rates[14] += fwd_rates[1154];

  //rxn 1155
  //sp 49
  sp_rates[48] -= fwd_rates[1155];
  //sp 14
  sp_rates[13] -= fwd_rates[1155];
  //sp 15
  sp_rates[14] += fwd_rates[1155];
  //sp 48
  sp_rates[47] += fwd_rates[1155];

  //rxn 1156
  //sp 36
  sp_rates[35] -= fwd_rates[1156];
  //sp 37
  sp_rates[36] += fwd_rates[1156];
  //sp 14
  sp_rates[13] += fwd_rates[1156];
  //sp 15
  sp_rates[14] -= fwd_rates[1156];

  //rxn 1157
  //sp 34
  sp_rates[33] += fwd_rates[1157];
  //sp 35
  sp_rates[34] -= fwd_rates[1157];
  //sp 36
  sp_rates[35] -= fwd_rates[1157];
  //sp 37
  sp_rates[36] += fwd_rates[1157];

  //rxn 1158
  //sp 36
  sp_rates[35] -= fwd_rates[1158];
  //sp 37
  sp_rates[36] += fwd_rates[1158];
  //sp 30
  sp_rates[29] += fwd_rates[1158];
  //sp 31
  sp_rates[30] -= fwd_rates[1158];

  //rxn 1159
  //sp 59
  sp_rates[58] += fwd_rates[1159];
  //sp 36
  sp_rates[35] -= fwd_rates[1159];
  //sp 61
  sp_rates[60] -= fwd_rates[1159];
  //sp 37
  sp_rates[36] += fwd_rates[1159];

  //rxn 1160
  //sp 60
  sp_rates[59] -= fwd_rates[1160];
  //sp 59
  sp_rates[58] += fwd_rates[1160];
  //sp 36
  sp_rates[35] -= fwd_rates[1160];
  //sp 37
  sp_rates[36] += fwd_rates[1160];

  //rxn 1161
  //sp 36
  sp_rates[35] -= fwd_rates[1161];
  //sp 37
  sp_rates[36] += fwd_rates[1161];
  //sp 62
  sp_rates[61] += fwd_rates[1161];
  //sp 64
  sp_rates[63] -= fwd_rates[1161];

  //rxn 1162
  //sp 65
  sp_rates[64] -= fwd_rates[1162];
  //sp 36
  sp_rates[35] -= fwd_rates[1162];
  //sp 37
  sp_rates[36] += fwd_rates[1162];
  //sp 62
  sp_rates[61] += fwd_rates[1162];

  //rxn 1163
  //sp 36
  sp_rates[35] -= fwd_rates[1163];
  //sp 37
  sp_rates[36] += fwd_rates[1163];
  //sp 62
  sp_rates[61] += fwd_rates[1163];
  //sp 63
  sp_rates[62] -= fwd_rates[1163];

  //rxn 1164
  //sp 36
  sp_rates[35] -= fwd_rates[1164];
  //sp 37
  sp_rates[36] += fwd_rates[1164];
  //sp 102
  sp_rates[101] += fwd_rates[1164];
  //sp 103
  sp_rates[102] -= fwd_rates[1164];

  //rxn 1165
  //sp 100
  sp_rates[99] -= fwd_rates[1165];
  //sp 36
  sp_rates[35] -= fwd_rates[1165];
  //sp 37
  sp_rates[36] += fwd_rates[1165];
  //sp 97
  sp_rates[96] += fwd_rates[1165];

  //rxn 1166
  //sp 98
  sp_rates[97] += fwd_rates[1166];
  //sp 99
  sp_rates[98] -= fwd_rates[1166];
  //sp 36
  sp_rates[35] -= fwd_rates[1166];
  //sp 37
  sp_rates[36] += fwd_rates[1166];

  //rxn 1167
  //sp 36
  sp_rates[35] -= fwd_rates[1167];
  //sp 37
  sp_rates[36] += fwd_rates[1167];
  //sp 94
  sp_rates[93] += fwd_rates[1167];
  //sp 95
  sp_rates[94] -= fwd_rates[1167];

  //rxn 1168
  //sp 105
  sp_rates[104] -= fwd_rates[1168];
  //sp 36
  sp_rates[35] -= fwd_rates[1168];
  //sp 37
  sp_rates[36] += fwd_rates[1168];
  //sp 104
  sp_rates[103] += fwd_rates[1168];

  //rxn 1169
  //sp 107
  sp_rates[106] += fwd_rates[1169];
  //sp 36
  sp_rates[35] -= fwd_rates[1169];
  //sp 109
  sp_rates[108] -= fwd_rates[1169];
  //sp 37
  sp_rates[36] += fwd_rates[1169];

  //rxn 1170
  //sp 36
  sp_rates[35] -= fwd_rates[1170];
  //sp 37
  sp_rates[36] += fwd_rates[1170];
  //sp 110
  sp_rates[109] += fwd_rates[1170];
  //sp 111
  sp_rates[110] -= fwd_rates[1170];

  //rxn 1171
  //sp 36
  sp_rates[35] -= fwd_rates[1171];
  //sp 37
  sp_rates[36] += fwd_rates[1171];
  //sp 6
  sp_rates[5] -= fwd_rates[1171];
  //sp 8
  sp_rates[7] += fwd_rates[1171];

  //rxn 1172
  //sp 10
  sp_rates[9] += fwd_rates[1172];
  //sp 36
  sp_rates[35] -= fwd_rates[1172];
  //sp 5
  sp_rates[4] -= fwd_rates[1172];
  //sp 37
  sp_rates[36] += fwd_rates[1172];

  //rxn 1173
  //sp 9
  sp_rates[8] += fwd_rates[1173];
  //sp 10
  sp_rates[9] -= fwd_rates[1173];
  //sp 36
  sp_rates[35] -= fwd_rates[1173];
  //sp 37
  sp_rates[36] += fwd_rates[1173];

  //rxn 1174
  //sp 36
  sp_rates[35] -= fwd_rates[1174];
  //sp 37
  sp_rates[36] += fwd_rates[1174];
  //sp 22
  sp_rates[21] += fwd_rates[1174];
  //sp 23
  sp_rates[22] -= fwd_rates[1174];

  //rxn 1175
  //sp 36
  sp_rates[35] -= fwd_rates[1175];
  //sp 37
  sp_rates[36] += fwd_rates[1175];
  //sp 22
  sp_rates[21] += fwd_rates[1175];
  //sp 24
  sp_rates[23] -= fwd_rates[1175];

  //rxn 1176
  //sp 25
  sp_rates[24] += fwd_rates[1176];
  //sp 26
  sp_rates[25] -= fwd_rates[1176];
  //sp 36
  sp_rates[35] -= fwd_rates[1176];
  //sp 37
  sp_rates[36] += fwd_rates[1176];

  //rxn 1177
  //sp 36
  sp_rates[35] -= fwd_rates[1177];
  //sp 37
  sp_rates[36] += fwd_rates[1177];
  //sp 46
  sp_rates[45] -= fwd_rates[1177];
  //sp 45
  sp_rates[44] += fwd_rates[1177];

  //rxn 1178
  //sp 36
  sp_rates[35] -= fwd_rates[1178];
  //sp 37
  sp_rates[36] += fwd_rates[1178];
  //sp 47
  sp_rates[46] -= fwd_rates[1178];
  //sp 45
  sp_rates[44] += fwd_rates[1178];

  //rxn 1179
  //sp 49
  sp_rates[48] -= fwd_rates[1179];
  //sp 36
  sp_rates[35] -= fwd_rates[1179];
  //sp 37
  sp_rates[36] += fwd_rates[1179];
  //sp 48
  sp_rates[47] += fwd_rates[1179];

  //rxn 1180
  //sp 34
  sp_rates[33] += fwd_rates[1180];
  //sp 35
  sp_rates[34] -= fwd_rates[1180];
  //sp 30
  sp_rates[29] -= fwd_rates[1180];
  //sp 31
  sp_rates[30] += fwd_rates[1180];

  //rxn 1181
  //sp 62
  sp_rates[61] += fwd_rates[1181];
  //sp 30
  sp_rates[29] -= fwd_rates[1181];
  //sp 31
  sp_rates[30] += fwd_rates[1181];
  //sp 64
  sp_rates[63] -= fwd_rates[1181];

  //rxn 1182
  //sp 65
  sp_rates[64] -= fwd_rates[1182];
  //sp 62
  sp_rates[61] += fwd_rates[1182];
  //sp 30
  sp_rates[29] -= fwd_rates[1182];
  //sp 31
  sp_rates[30] += fwd_rates[1182];

  //rxn 1183
  //sp 102
  sp_rates[101] += fwd_rates[1183];
  //sp 30
  sp_rates[29] -= fwd_rates[1183];
  //sp 103
  sp_rates[102] -= fwd_rates[1183];
  //sp 31
  sp_rates[30] += fwd_rates[1183];

  //rxn 1184
  //sp 97
  sp_rates[96] += fwd_rates[1184];
  //sp 100
  sp_rates[99] -= fwd_rates[1184];
  //sp 30
  sp_rates[29] -= fwd_rates[1184];
  //sp 31
  sp_rates[30] += fwd_rates[1184];

  //rxn 1185
  //sp 98
  sp_rates[97] += fwd_rates[1185];
  //sp 99
  sp_rates[98] -= fwd_rates[1185];
  //sp 30
  sp_rates[29] -= fwd_rates[1185];
  //sp 31
  sp_rates[30] += fwd_rates[1185];

  //rxn 1186
  //sp 94
  sp_rates[93] += fwd_rates[1186];
  //sp 30
  sp_rates[29] -= fwd_rates[1186];
  //sp 95
  sp_rates[94] -= fwd_rates[1186];
  //sp 31
  sp_rates[30] += fwd_rates[1186];

  //rxn 1187
  //sp 105
  sp_rates[104] -= fwd_rates[1187];
  //sp 30
  sp_rates[29] -= fwd_rates[1187];
  //sp 31
  sp_rates[30] += fwd_rates[1187];
  //sp 104
  sp_rates[103] += fwd_rates[1187];

  //rxn 1188
  //sp 107
  sp_rates[106] += fwd_rates[1188];
  //sp 109
  sp_rates[108] -= fwd_rates[1188];
  //sp 30
  sp_rates[29] -= fwd_rates[1188];
  //sp 31
  sp_rates[30] += fwd_rates[1188];

  //rxn 1189
  //sp 110
  sp_rates[109] += fwd_rates[1189];
  //sp 30
  sp_rates[29] -= fwd_rates[1189];
  //sp 111
  sp_rates[110] -= fwd_rates[1189];
  //sp 31
  sp_rates[30] += fwd_rates[1189];

  //rxn 1190
  //sp 22
  sp_rates[21] += fwd_rates[1190];
  //sp 30
  sp_rates[29] -= fwd_rates[1190];
  //sp 31
  sp_rates[30] += fwd_rates[1190];
  //sp 24
  sp_rates[23] -= fwd_rates[1190];

  //rxn 1191
  //sp 25
  sp_rates[24] += fwd_rates[1191];
  //sp 26
  sp_rates[25] -= fwd_rates[1191];
  //sp 30
  sp_rates[29] -= fwd_rates[1191];
  //sp 31
  sp_rates[30] += fwd_rates[1191];

  //rxn 1192
  //sp 46
  sp_rates[45] -= fwd_rates[1192];
  //sp 45
  sp_rates[44] += fwd_rates[1192];
  //sp 30
  sp_rates[29] -= fwd_rates[1192];
  //sp 31
  sp_rates[30] += fwd_rates[1192];

  //rxn 1193
  //sp 45
  sp_rates[44] += fwd_rates[1193];
  //sp 30
  sp_rates[29] -= fwd_rates[1193];
  //sp 47
  sp_rates[46] -= fwd_rates[1193];
  //sp 31
  sp_rates[30] += fwd_rates[1193];

  //rxn 1194
  //sp 79
  sp_rates[78] += fwd_rates[1194];
  //sp 30
  sp_rates[29] -= fwd_rates[1194];
  //sp 31
  sp_rates[30] += fwd_rates[1194];
  //sp 80
  sp_rates[79] -= fwd_rates[1194];

  //rxn 1195
  //sp 49
  sp_rates[48] -= fwd_rates[1195];
  //sp 30
  sp_rates[29] -= fwd_rates[1195];
  //sp 31
  sp_rates[30] += fwd_rates[1195];
  //sp 48
  sp_rates[47] += fwd_rates[1195];

  //rxn 1196
  //sp 34
  sp_rates[33] -= fwd_rates[1196];
  //sp 35
  sp_rates[34] += fwd_rates[1196];
  //sp 30
  sp_rates[29] += fwd_rates[1196];
  //sp 31
  sp_rates[30] -= fwd_rates[1196];

  //rxn 1197
  //sp 34
  sp_rates[33] -= fwd_rates[1197];
  //sp 35
  sp_rates[34] += fwd_rates[1197];
  //sp 62
  sp_rates[61] += fwd_rates[1197];
  //sp 64
  sp_rates[63] -= fwd_rates[1197];

  //rxn 1198
  //sp 65
  sp_rates[64] -= fwd_rates[1198];
  //sp 34
  sp_rates[33] -= fwd_rates[1198];
  //sp 35
  sp_rates[34] += fwd_rates[1198];
  //sp 62
  sp_rates[61] += fwd_rates[1198];

  //rxn 1199
  //sp 34
  sp_rates[33] -= fwd_rates[1199];
  //sp 35
  sp_rates[34] += fwd_rates[1199];
  //sp 62
  sp_rates[61] += fwd_rates[1199];
  //sp 63
  sp_rates[62] -= fwd_rates[1199];

  //rxn 1200
  //sp 97
  sp_rates[96] += fwd_rates[1200];
  //sp 34
  sp_rates[33] -= fwd_rates[1200];
  //sp 35
  sp_rates[34] += fwd_rates[1200];
  //sp 100
  sp_rates[99] -= fwd_rates[1200];

  //rxn 1201
  //sp 34
  sp_rates[33] -= fwd_rates[1201];
  //sp 99
  sp_rates[98] -= fwd_rates[1201];
  //sp 35
  sp_rates[34] += fwd_rates[1201];
  //sp 98
  sp_rates[97] += fwd_rates[1201];

  //rxn 1202
  //sp 34
  sp_rates[33] -= fwd_rates[1202];
  //sp 35
  sp_rates[34] += fwd_rates[1202];
  //sp 94
  sp_rates[93] += fwd_rates[1202];
  //sp 95
  sp_rates[94] -= fwd_rates[1202];

  //rxn 1203
  //sp 105
  sp_rates[104] -= fwd_rates[1203];
  //sp 34
  sp_rates[33] -= fwd_rates[1203];
  //sp 35
  sp_rates[34] += fwd_rates[1203];
  //sp 104
  sp_rates[103] += fwd_rates[1203];

  //rxn 1204
  //sp 34
  sp_rates[33] -= fwd_rates[1204];
  //sp 35
  sp_rates[34] += fwd_rates[1204];
  //sp 109
  sp_rates[108] -= fwd_rates[1204];
  //sp 107
  sp_rates[106] += fwd_rates[1204];

  //rxn 1205
  //sp 34
  sp_rates[33] -= fwd_rates[1205];
  //sp 35
  sp_rates[34] += fwd_rates[1205];
  //sp 6
  sp_rates[5] -= fwd_rates[1205];
  //sp 8
  sp_rates[7] += fwd_rates[1205];

  //rxn 1206
  //sp 9
  sp_rates[8] += fwd_rates[1206];
  //sp 10
  sp_rates[9] -= fwd_rates[1206];
  //sp 35
  sp_rates[34] += fwd_rates[1206];
  //sp 34
  sp_rates[33] -= fwd_rates[1206];

  //rxn 1207
  //sp 34
  sp_rates[33] -= fwd_rates[1207];
  //sp 35
  sp_rates[34] += fwd_rates[1207];
  //sp 22
  sp_rates[21] += fwd_rates[1207];
  //sp 24
  sp_rates[23] -= fwd_rates[1207];

  //rxn 1208
  //sp 34
  sp_rates[33] -= fwd_rates[1208];
  //sp 35
  sp_rates[34] += fwd_rates[1208];
  //sp 45
  sp_rates[44] += fwd_rates[1208];
  //sp 46
  sp_rates[45] -= fwd_rates[1208];

  //rxn 1209
  //sp 34
  sp_rates[33] -= fwd_rates[1209];
  //sp 35
  sp_rates[34] += fwd_rates[1209];
  //sp 45
  sp_rates[44] += fwd_rates[1209];
  //sp 47
  sp_rates[46] -= fwd_rates[1209];

  //rxn 1210
  //sp 34
  sp_rates[33] -= fwd_rates[1210];
  //sp 35
  sp_rates[34] += fwd_rates[1210];
  //sp 79
  sp_rates[78] += fwd_rates[1210];
  //sp 80
  sp_rates[79] -= fwd_rates[1210];

  //rxn 1211
  //sp 49
  sp_rates[48] -= fwd_rates[1211];
  //sp 34
  sp_rates[33] -= fwd_rates[1211];
  //sp 35
  sp_rates[34] += fwd_rates[1211];
  //sp 48
  sp_rates[47] += fwd_rates[1211];

  //rxn 1212
  //sp 59
  sp_rates[58] -= fwd_rates[1212];
  //sp 61
  sp_rates[60] += fwd_rates[1212];
  //sp 62
  sp_rates[61] += fwd_rates[1212];
  //sp 64
  sp_rates[63] -= fwd_rates[1212];

  //rxn 1213
  //sp 65
  sp_rates[64] -= fwd_rates[1213];
  //sp 59
  sp_rates[58] -= fwd_rates[1213];
  //sp 61
  sp_rates[60] += fwd_rates[1213];
  //sp 62
  sp_rates[61] += fwd_rates[1213];

  //rxn 1214
  //sp 59
  sp_rates[58] -= fwd_rates[1214];
  //sp 61
  sp_rates[60] += fwd_rates[1214];
  //sp 102
  sp_rates[101] += fwd_rates[1214];
  //sp 103
  sp_rates[102] -= fwd_rates[1214];

  //rxn 1215
  //sp 97
  sp_rates[96] += fwd_rates[1215];
  //sp 59
  sp_rates[58] -= fwd_rates[1215];
  //sp 100
  sp_rates[99] -= fwd_rates[1215];
  //sp 61
  sp_rates[60] += fwd_rates[1215];

  //rxn 1216
  //sp 98
  sp_rates[97] += fwd_rates[1216];
  //sp 59
  sp_rates[58] -= fwd_rates[1216];
  //sp 61
  sp_rates[60] += fwd_rates[1216];
  //sp 99
  sp_rates[98] -= fwd_rates[1216];

  //rxn 1217
  //sp 59
  sp_rates[58] -= fwd_rates[1217];
  //sp 61
  sp_rates[60] += fwd_rates[1217];
  //sp 94
  sp_rates[93] += fwd_rates[1217];
  //sp 95
  sp_rates[94] -= fwd_rates[1217];

  //rxn 1218
  //sp 105
  sp_rates[104] -= fwd_rates[1218];
  //sp 59
  sp_rates[58] -= fwd_rates[1218];
  //sp 61
  sp_rates[60] += fwd_rates[1218];
  //sp 104
  sp_rates[103] += fwd_rates[1218];

  //rxn 1219
  //sp 59
  sp_rates[58] -= fwd_rates[1219];
  //sp 109
  sp_rates[108] -= fwd_rates[1219];
  //sp 61
  sp_rates[60] += fwd_rates[1219];
  //sp 107
  sp_rates[106] += fwd_rates[1219];

  //rxn 1220
  //sp 59
  sp_rates[58] -= fwd_rates[1220];
  //sp 61
  sp_rates[60] += fwd_rates[1220];
  //sp 110
  sp_rates[109] += fwd_rates[1220];
  //sp 111
  sp_rates[110] -= fwd_rates[1220];

  //rxn 1221
  //sp 59
  sp_rates[58] -= fwd_rates[1221];
  //sp 61
  sp_rates[60] += fwd_rates[1221];
  //sp 22
  sp_rates[21] += fwd_rates[1221];
  //sp 24
  sp_rates[23] -= fwd_rates[1221];

  //rxn 1222
  //sp 25
  sp_rates[24] += fwd_rates[1222];
  //sp 26
  sp_rates[25] -= fwd_rates[1222];
  //sp 59
  sp_rates[58] -= fwd_rates[1222];
  //sp 61
  sp_rates[60] += fwd_rates[1222];

  //rxn 1223
  //sp 59
  sp_rates[58] -= fwd_rates[1223];
  //sp 45
  sp_rates[44] += fwd_rates[1223];
  //sp 46
  sp_rates[45] -= fwd_rates[1223];
  //sp 61
  sp_rates[60] += fwd_rates[1223];

  //rxn 1224
  //sp 59
  sp_rates[58] -= fwd_rates[1224];
  //sp 45
  sp_rates[44] += fwd_rates[1224];
  //sp 47
  sp_rates[46] -= fwd_rates[1224];
  //sp 61
  sp_rates[60] += fwd_rates[1224];

  //rxn 1225
  //sp 59
  sp_rates[58] -= fwd_rates[1225];
  //sp 61
  sp_rates[60] += fwd_rates[1225];
  //sp 79
  sp_rates[78] += fwd_rates[1225];
  //sp 80
  sp_rates[79] -= fwd_rates[1225];

  //rxn 1226
  //sp 49
  sp_rates[48] -= fwd_rates[1226];
  //sp 59
  sp_rates[58] -= fwd_rates[1226];
  //sp 61
  sp_rates[60] += fwd_rates[1226];
  //sp 48
  sp_rates[47] += fwd_rates[1226];

  //rxn 1227
  //sp 59
  sp_rates[58] -= fwd_rates[1227];
  //sp 60
  sp_rates[59] += fwd_rates[1227];
  //sp 62
  sp_rates[61] += fwd_rates[1227];
  //sp 64
  sp_rates[63] -= fwd_rates[1227];

  //rxn 1228
  //sp 65
  sp_rates[64] -= fwd_rates[1228];
  //sp 59
  sp_rates[58] -= fwd_rates[1228];
  //sp 60
  sp_rates[59] += fwd_rates[1228];
  //sp 62
  sp_rates[61] += fwd_rates[1228];

  //rxn 1229
  //sp 59
  sp_rates[58] -= fwd_rates[1229];
  //sp 60
  sp_rates[59] += fwd_rates[1229];
  //sp 102
  sp_rates[101] += fwd_rates[1229];
  //sp 103
  sp_rates[102] -= fwd_rates[1229];

  //rxn 1230
  //sp 97
  sp_rates[96] += fwd_rates[1230];
  //sp 60
  sp_rates[59] += fwd_rates[1230];
  //sp 59
  sp_rates[58] -= fwd_rates[1230];
  //sp 100
  sp_rates[99] -= fwd_rates[1230];

  //rxn 1231
  //sp 98
  sp_rates[97] += fwd_rates[1231];
  //sp 59
  sp_rates[58] -= fwd_rates[1231];
  //sp 60
  sp_rates[59] += fwd_rates[1231];
  //sp 99
  sp_rates[98] -= fwd_rates[1231];

  //rxn 1232
  //sp 59
  sp_rates[58] -= fwd_rates[1232];
  //sp 60
  sp_rates[59] += fwd_rates[1232];
  //sp 94
  sp_rates[93] += fwd_rates[1232];
  //sp 95
  sp_rates[94] -= fwd_rates[1232];

  //rxn 1233
  //sp 105
  sp_rates[104] -= fwd_rates[1233];
  //sp 59
  sp_rates[58] -= fwd_rates[1233];
  //sp 60
  sp_rates[59] += fwd_rates[1233];
  //sp 104
  sp_rates[103] += fwd_rates[1233];

  //rxn 1234
  //sp 59
  sp_rates[58] -= fwd_rates[1234];
  //sp 60
  sp_rates[59] += fwd_rates[1234];
  //sp 109
  sp_rates[108] -= fwd_rates[1234];
  //sp 107
  sp_rates[106] += fwd_rates[1234];

  //rxn 1235
  //sp 59
  sp_rates[58] -= fwd_rates[1235];
  //sp 60
  sp_rates[59] += fwd_rates[1235];
  //sp 110
  sp_rates[109] += fwd_rates[1235];
  //sp 111
  sp_rates[110] -= fwd_rates[1235];

  //rxn 1236
  //sp 59
  sp_rates[58] -= fwd_rates[1236];
  //sp 60
  sp_rates[59] += fwd_rates[1236];
  //sp 22
  sp_rates[21] += fwd_rates[1236];
  //sp 24
  sp_rates[23] -= fwd_rates[1236];

  //rxn 1237
  //sp 25
  sp_rates[24] += fwd_rates[1237];
  //sp 26
  sp_rates[25] -= fwd_rates[1237];
  //sp 59
  sp_rates[58] -= fwd_rates[1237];
  //sp 60
  sp_rates[59] += fwd_rates[1237];

  //rxn 1238
  //sp 59
  sp_rates[58] -= fwd_rates[1238];
  //sp 60
  sp_rates[59] += fwd_rates[1238];
  //sp 45
  sp_rates[44] += fwd_rates[1238];
  //sp 46
  sp_rates[45] -= fwd_rates[1238];

  //rxn 1239
  //sp 59
  sp_rates[58] -= fwd_rates[1239];
  //sp 60
  sp_rates[59] += fwd_rates[1239];
  //sp 45
  sp_rates[44] += fwd_rates[1239];
  //sp 47
  sp_rates[46] -= fwd_rates[1239];

  //rxn 1240
  //sp 59
  sp_rates[58] -= fwd_rates[1240];
  //sp 60
  sp_rates[59] += fwd_rates[1240];
  //sp 79
  sp_rates[78] += fwd_rates[1240];
  //sp 80
  sp_rates[79] -= fwd_rates[1240];

  //rxn 1241
  //sp 49
  sp_rates[48] -= fwd_rates[1241];
  //sp 59
  sp_rates[58] -= fwd_rates[1241];
  //sp 60
  sp_rates[59] += fwd_rates[1241];
  //sp 48
  sp_rates[47] += fwd_rates[1241];

  //rxn 1242
  //sp 34
  sp_rates[33] += fwd_rates[1242];
  //sp 35
  sp_rates[34] -= fwd_rates[1242];
  //sp 62
  sp_rates[61] -= fwd_rates[1242];
  //sp 63
  sp_rates[62] += fwd_rates[1242];

  //rxn 1243
  //sp 63
  sp_rates[62] += fwd_rates[1243];
  //sp 64
  sp_rates[63] -= fwd_rates[1243];

  //rxn 1244
  //sp 65
  sp_rates[64] -= fwd_rates[1244];
  //sp 63
  sp_rates[62] += fwd_rates[1244];

  //rxn 1245
  //sp 97
  sp_rates[96] += fwd_rates[1245];
  //sp 100
  sp_rates[99] -= fwd_rates[1245];
  //sp 62
  sp_rates[61] -= fwd_rates[1245];
  //sp 63
  sp_rates[62] += fwd_rates[1245];

  //rxn 1246
  //sp 98
  sp_rates[97] += fwd_rates[1246];
  //sp 99
  sp_rates[98] -= fwd_rates[1246];
  //sp 62
  sp_rates[61] -= fwd_rates[1246];
  //sp 63
  sp_rates[62] += fwd_rates[1246];

  //rxn 1247
  //sp 105
  sp_rates[104] -= fwd_rates[1247];
  //sp 62
  sp_rates[61] -= fwd_rates[1247];
  //sp 63
  sp_rates[62] += fwd_rates[1247];
  //sp 104
  sp_rates[103] += fwd_rates[1247];

  //rxn 1248
  //sp 107
  sp_rates[106] += fwd_rates[1248];
  //sp 109
  sp_rates[108] -= fwd_rates[1248];
  //sp 62
  sp_rates[61] -= fwd_rates[1248];
  //sp 63
  sp_rates[62] += fwd_rates[1248];

  //rxn 1249
  //sp 63
  sp_rates[62] += fwd_rates[1249];
  //sp 110
  sp_rates[109] += fwd_rates[1249];
  //sp 62
  sp_rates[61] -= fwd_rates[1249];
  //sp 111
  sp_rates[110] -= fwd_rates[1249];

  //rxn 1250
  //sp 22
  sp_rates[21] += fwd_rates[1250];
  //sp 62
  sp_rates[61] -= fwd_rates[1250];
  //sp 63
  sp_rates[62] += fwd_rates[1250];
  //sp 24
  sp_rates[23] -= fwd_rates[1250];

  //rxn 1251
  //sp 25
  sp_rates[24] += fwd_rates[1251];
  //sp 26
  sp_rates[25] -= fwd_rates[1251];
  //sp 62
  sp_rates[61] -= fwd_rates[1251];
  //sp 63
  sp_rates[62] += fwd_rates[1251];

  //rxn 1252
  //sp 62
  sp_rates[61] -= fwd_rates[1252];
  //sp 45
  sp_rates[44] += fwd_rates[1252];
  //sp 46
  sp_rates[45] -= fwd_rates[1252];
  //sp 63
  sp_rates[62] += fwd_rates[1252];

  //rxn 1253
  //sp 63
  sp_rates[62] += fwd_rates[1253];
  //sp 45
  sp_rates[44] += fwd_rates[1253];
  //sp 62
  sp_rates[61] -= fwd_rates[1253];
  //sp 47
  sp_rates[46] -= fwd_rates[1253];

  //rxn 1254
  //sp 49
  sp_rates[48] -= fwd_rates[1254];
  //sp 62
  sp_rates[61] -= fwd_rates[1254];
  //sp 63
  sp_rates[62] += fwd_rates[1254];
  //sp 48
  sp_rates[47] += fwd_rates[1254];

  //rxn 1255
  //sp 34
  sp_rates[33] += fwd_rates[1255];
  //sp 35
  sp_rates[34] -= fwd_rates[1255];
  //sp 62
  sp_rates[61] -= fwd_rates[1255];
  //sp 64
  sp_rates[63] += fwd_rates[1255];

  //rxn 1256
  //sp 30
  sp_rates[29] += fwd_rates[1256];
  //sp 62
  sp_rates[61] -= fwd_rates[1256];
  //sp 31
  sp_rates[30] -= fwd_rates[1256];
  //sp 64
  sp_rates[63] += fwd_rates[1256];

  //rxn 1257
  //sp 59
  sp_rates[58] += fwd_rates[1257];
  //sp 61
  sp_rates[60] -= fwd_rates[1257];
  //sp 62
  sp_rates[61] -= fwd_rates[1257];
  //sp 64
  sp_rates[63] += fwd_rates[1257];

  //rxn 1258
  //sp 59
  sp_rates[58] += fwd_rates[1258];
  //sp 60
  sp_rates[59] -= fwd_rates[1258];
  //sp 62
  sp_rates[61] -= fwd_rates[1258];
  //sp 64
  sp_rates[63] += fwd_rates[1258];

  //rxn 1259
  //sp 65
  sp_rates[64] -= fwd_rates[1259];
  //sp 64
  sp_rates[63] += fwd_rates[1259];

  //rxn 1260
  //sp 63
  sp_rates[62] -= fwd_rates[1260];
  //sp 64
  sp_rates[63] += fwd_rates[1260];

  //rxn 1261
  //sp 102
  sp_rates[101] += fwd_rates[1261];
  //sp 62
  sp_rates[61] -= fwd_rates[1261];
  //sp 103
  sp_rates[102] -= fwd_rates[1261];
  //sp 64
  sp_rates[63] += fwd_rates[1261];

  //rxn 1262
  //sp 97
  sp_rates[96] += fwd_rates[1262];
  //sp 100
  sp_rates[99] -= fwd_rates[1262];
  //sp 62
  sp_rates[61] -= fwd_rates[1262];
  //sp 64
  sp_rates[63] += fwd_rates[1262];

  //rxn 1263
  //sp 98
  sp_rates[97] += fwd_rates[1263];
  //sp 99
  sp_rates[98] -= fwd_rates[1263];
  //sp 62
  sp_rates[61] -= fwd_rates[1263];
  //sp 64
  sp_rates[63] += fwd_rates[1263];

  //rxn 1264
  //sp 105
  sp_rates[104] -= fwd_rates[1264];
  //sp 104
  sp_rates[103] += fwd_rates[1264];
  //sp 62
  sp_rates[61] -= fwd_rates[1264];
  //sp 64
  sp_rates[63] += fwd_rates[1264];

  //rxn 1265
  //sp 107
  sp_rates[106] += fwd_rates[1265];
  //sp 109
  sp_rates[108] -= fwd_rates[1265];
  //sp 62
  sp_rates[61] -= fwd_rates[1265];
  //sp 64
  sp_rates[63] += fwd_rates[1265];

  //rxn 1266
  //sp 110
  sp_rates[109] += fwd_rates[1266];
  //sp 62
  sp_rates[61] -= fwd_rates[1266];
  //sp 111
  sp_rates[110] -= fwd_rates[1266];
  //sp 64
  sp_rates[63] += fwd_rates[1266];

  //rxn 1267
  //sp 22
  sp_rates[21] += fwd_rates[1267];
  //sp 62
  sp_rates[61] -= fwd_rates[1267];
  //sp 23
  sp_rates[22] -= fwd_rates[1267];
  //sp 64
  sp_rates[63] += fwd_rates[1267];

  //rxn 1268
  //sp 22
  sp_rates[21] += fwd_rates[1268];
  //sp 64
  sp_rates[63] += fwd_rates[1268];
  //sp 62
  sp_rates[61] -= fwd_rates[1268];
  //sp 24
  sp_rates[23] -= fwd_rates[1268];

  //rxn 1269
  //sp 62
  sp_rates[61] -= fwd_rates[1269];
  //sp 45
  sp_rates[44] += fwd_rates[1269];
  //sp 46
  sp_rates[45] -= fwd_rates[1269];
  //sp 64
  sp_rates[63] += fwd_rates[1269];

  //rxn 1270
  //sp 45
  sp_rates[44] += fwd_rates[1270];
  //sp 62
  sp_rates[61] -= fwd_rates[1270];
  //sp 47
  sp_rates[46] -= fwd_rates[1270];
  //sp 64
  sp_rates[63] += fwd_rates[1270];

  //rxn 1271
  //sp 64
  sp_rates[63] += fwd_rates[1271];
  //sp 62
  sp_rates[61] -= fwd_rates[1271];
  //sp 79
  sp_rates[78] += fwd_rates[1271];
  //sp 80
  sp_rates[79] -= fwd_rates[1271];

  //rxn 1272
  //sp 49
  sp_rates[48] -= fwd_rates[1272];
  //sp 64
  sp_rates[63] += fwd_rates[1272];
  //sp 62
  sp_rates[61] -= fwd_rates[1272];
  //sp 48
  sp_rates[47] += fwd_rates[1272];

  //rxn 1273
  //sp 65
  sp_rates[64] += fwd_rates[1273];
  //sp 34
  sp_rates[33] += fwd_rates[1273];
  //sp 35
  sp_rates[34] -= fwd_rates[1273];
  //sp 62
  sp_rates[61] -= fwd_rates[1273];

  //rxn 1274
  //sp 65
  sp_rates[64] += fwd_rates[1274];
  //sp 30
  sp_rates[29] += fwd_rates[1274];
  //sp 62
  sp_rates[61] -= fwd_rates[1274];
  //sp 31
  sp_rates[30] -= fwd_rates[1274];

  //rxn 1275
  //sp 65
  sp_rates[64] += fwd_rates[1275];
  //sp 59
  sp_rates[58] += fwd_rates[1275];
  //sp 61
  sp_rates[60] -= fwd_rates[1275];
  //sp 62
  sp_rates[61] -= fwd_rates[1275];

  //rxn 1276
  //sp 65
  sp_rates[64] += fwd_rates[1276];
  //sp 59
  sp_rates[58] += fwd_rates[1276];
  //sp 60
  sp_rates[59] -= fwd_rates[1276];
  //sp 62
  sp_rates[61] -= fwd_rates[1276];

  //rxn 1277
  //sp 65
  sp_rates[64] += fwd_rates[1277];
  //sp 64
  sp_rates[63] -= fwd_rates[1277];

  //rxn 1278
  //sp 65
  sp_rates[64] += fwd_rates[1278];
  //sp 63
  sp_rates[62] -= fwd_rates[1278];

  //rxn 1279
  //sp 65
  sp_rates[64] += fwd_rates[1279];
  //sp 102
  sp_rates[101] += fwd_rates[1279];
  //sp 62
  sp_rates[61] -= fwd_rates[1279];
  //sp 103
  sp_rates[102] -= fwd_rates[1279];

  //rxn 1280
  //sp 65
  sp_rates[64] += fwd_rates[1280];
  //sp 100
  sp_rates[99] -= fwd_rates[1280];
  //sp 97
  sp_rates[96] += fwd_rates[1280];
  //sp 62
  sp_rates[61] -= fwd_rates[1280];

  //rxn 1281
  //sp 65
  sp_rates[64] += fwd_rates[1281];
  //sp 98
  sp_rates[97] += fwd_rates[1281];
  //sp 99
  sp_rates[98] -= fwd_rates[1281];
  //sp 62
  sp_rates[61] -= fwd_rates[1281];

  //rxn 1282
  //sp 105
  sp_rates[104] -= fwd_rates[1282];
  //sp 65
  sp_rates[64] += fwd_rates[1282];
  //sp 62
  sp_rates[61] -= fwd_rates[1282];
  //sp 104
  sp_rates[103] += fwd_rates[1282];

  //rxn 1283
  //sp 65
  sp_rates[64] += fwd_rates[1283];
  //sp 107
  sp_rates[106] += fwd_rates[1283];
  //sp 109
  sp_rates[108] -= fwd_rates[1283];
  //sp 62
  sp_rates[61] -= fwd_rates[1283];

  //rxn 1284
  //sp 65
  sp_rates[64] += fwd_rates[1284];
  //sp 110
  sp_rates[109] += fwd_rates[1284];
  //sp 62
  sp_rates[61] -= fwd_rates[1284];
  //sp 111
  sp_rates[110] -= fwd_rates[1284];

  //rxn 1285
  //sp 65
  sp_rates[64] += fwd_rates[1285];
  //sp 22
  sp_rates[21] += fwd_rates[1285];
  //sp 62
  sp_rates[61] -= fwd_rates[1285];
  //sp 23
  sp_rates[22] -= fwd_rates[1285];

  //rxn 1286
  //sp 65
  sp_rates[64] += fwd_rates[1286];
  //sp 22
  sp_rates[21] += fwd_rates[1286];
  //sp 62
  sp_rates[61] -= fwd_rates[1286];
  //sp 24
  sp_rates[23] -= fwd_rates[1286];

  //rxn 1287
  //sp 65
  sp_rates[64] += fwd_rates[1287];
  //sp 62
  sp_rates[61] -= fwd_rates[1287];
  //sp 45
  sp_rates[44] += fwd_rates[1287];
  //sp 46
  sp_rates[45] -= fwd_rates[1287];

  //rxn 1288
  //sp 65
  sp_rates[64] += fwd_rates[1288];
  //sp 45
  sp_rates[44] += fwd_rates[1288];
  //sp 62
  sp_rates[61] -= fwd_rates[1288];
  //sp 47
  sp_rates[46] -= fwd_rates[1288];

  //rxn 1289
  //sp 65
  sp_rates[64] += fwd_rates[1289];
  //sp 62
  sp_rates[61] -= fwd_rates[1289];
  //sp 79
  sp_rates[78] += fwd_rates[1289];
  //sp 80
  sp_rates[79] -= fwd_rates[1289];

  //rxn 1290
  //sp 49
  sp_rates[48] -= fwd_rates[1290];
  //sp 65
  sp_rates[64] += fwd_rates[1290];
  //sp 62
  sp_rates[61] -= fwd_rates[1290];
  //sp 48
  sp_rates[47] += fwd_rates[1290];

  //rxn 1291
  //sp 34
  sp_rates[33] += fwd_rates[1291];
  //sp 35
  sp_rates[34] -= fwd_rates[1291];
  //sp 79
  sp_rates[78] -= fwd_rates[1291];
  //sp 80
  sp_rates[79] += fwd_rates[1291];

  //rxn 1292
  //sp 79
  sp_rates[78] -= fwd_rates[1292];
  //sp 30
  sp_rates[29] += fwd_rates[1292];
  //sp 31
  sp_rates[30] -= fwd_rates[1292];
  //sp 80
  sp_rates[79] += fwd_rates[1292];

  //rxn 1293
  //sp 59
  sp_rates[58] += fwd_rates[1293];
  //sp 61
  sp_rates[60] -= fwd_rates[1293];
  //sp 79
  sp_rates[78] -= fwd_rates[1293];
  //sp 80
  sp_rates[79] += fwd_rates[1293];

  //rxn 1294
  //sp 59
  sp_rates[58] += fwd_rates[1294];
  //sp 60
  sp_rates[59] -= fwd_rates[1294];
  //sp 79
  sp_rates[78] -= fwd_rates[1294];
  //sp 80
  sp_rates[79] += fwd_rates[1294];

  //rxn 1295
  //sp 80
  sp_rates[79] += fwd_rates[1295];
  //sp 62
  sp_rates[61] += fwd_rates[1295];
  //sp 79
  sp_rates[78] -= fwd_rates[1295];
  //sp 64
  sp_rates[63] -= fwd_rates[1295];

  //rxn 1296
  //sp 65
  sp_rates[64] -= fwd_rates[1296];
  //sp 62
  sp_rates[61] += fwd_rates[1296];
  //sp 79
  sp_rates[78] -= fwd_rates[1296];
  //sp 80
  sp_rates[79] += fwd_rates[1296];

  //rxn 1297
  //sp 97
  sp_rates[96] += fwd_rates[1297];
  //sp 100
  sp_rates[99] -= fwd_rates[1297];
  //sp 79
  sp_rates[78] -= fwd_rates[1297];
  //sp 80
  sp_rates[79] += fwd_rates[1297];

  //rxn 1298
  //sp 98
  sp_rates[97] += fwd_rates[1298];
  //sp 99
  sp_rates[98] -= fwd_rates[1298];
  //sp 79
  sp_rates[78] -= fwd_rates[1298];
  //sp 80
  sp_rates[79] += fwd_rates[1298];

  //rxn 1299
  //sp 95
  sp_rates[94] -= fwd_rates[1299];
  //sp 94
  sp_rates[93] += fwd_rates[1299];
  //sp 79
  sp_rates[78] -= fwd_rates[1299];
  //sp 80
  sp_rates[79] += fwd_rates[1299];

  //rxn 1300
  //sp 105
  sp_rates[104] -= fwd_rates[1300];
  //sp 104
  sp_rates[103] += fwd_rates[1300];
  //sp 79
  sp_rates[78] -= fwd_rates[1300];
  //sp 80
  sp_rates[79] += fwd_rates[1300];

  //rxn 1301
  //sp 107
  sp_rates[106] += fwd_rates[1301];
  //sp 109
  sp_rates[108] -= fwd_rates[1301];
  //sp 79
  sp_rates[78] -= fwd_rates[1301];
  //sp 80
  sp_rates[79] += fwd_rates[1301];

  //rxn 1302
  //sp 111
  sp_rates[110] -= fwd_rates[1302];
  //sp 110
  sp_rates[109] += fwd_rates[1302];
  //sp 79
  sp_rates[78] -= fwd_rates[1302];
  //sp 80
  sp_rates[79] += fwd_rates[1302];

  //rxn 1303
  //sp 80
  sp_rates[79] += fwd_rates[1303];
  //sp 6
  sp_rates[5] -= fwd_rates[1303];
  //sp 79
  sp_rates[78] -= fwd_rates[1303];
  //sp 8
  sp_rates[7] += fwd_rates[1303];

  //rxn 1304
  //sp 79
  sp_rates[78] -= fwd_rates[1304];
  //sp 22
  sp_rates[21] += fwd_rates[1304];
  //sp 23
  sp_rates[22] -= fwd_rates[1304];
  //sp 80
  sp_rates[79] += fwd_rates[1304];

  //rxn 1305
  //sp 80
  sp_rates[79] += fwd_rates[1305];
  //sp 22
  sp_rates[21] += fwd_rates[1305];
  //sp 79
  sp_rates[78] -= fwd_rates[1305];
  //sp 24
  sp_rates[23] -= fwd_rates[1305];

  //rxn 1306
  //sp 25
  sp_rates[24] += fwd_rates[1306];
  //sp 26
  sp_rates[25] -= fwd_rates[1306];
  //sp 79
  sp_rates[78] -= fwd_rates[1306];
  //sp 80
  sp_rates[79] += fwd_rates[1306];

  //rxn 1307
  //sp 45
  sp_rates[44] += fwd_rates[1307];
  //sp 46
  sp_rates[45] -= fwd_rates[1307];
  //sp 79
  sp_rates[78] -= fwd_rates[1307];
  //sp 80
  sp_rates[79] += fwd_rates[1307];

  //rxn 1308
  //sp 79
  sp_rates[78] -= fwd_rates[1308];
  //sp 45
  sp_rates[44] += fwd_rates[1308];
  //sp 47
  sp_rates[46] -= fwd_rates[1308];
  //sp 80
  sp_rates[79] += fwd_rates[1308];

  //rxn 1309
  //sp 49
  sp_rates[48] -= fwd_rates[1309];
  //sp 80
  sp_rates[79] += fwd_rates[1309];
  //sp 79
  sp_rates[78] -= fwd_rates[1309];
  //sp 48
  sp_rates[47] += fwd_rates[1309];

  //rxn 1310
  //sp 30
  sp_rates[29] += fwd_rates[1310];
  //sp 78
  sp_rates[77] -= fwd_rates[1310];
  //sp 31
  sp_rates[30] -= fwd_rates[1310];
  //sp 80
  sp_rates[79] += fwd_rates[1310];

  //rxn 1311
  //sp 59
  sp_rates[58] += fwd_rates[1311];
  //sp 61
  sp_rates[60] -= fwd_rates[1311];
  //sp 78
  sp_rates[77] -= fwd_rates[1311];
  //sp 80
  sp_rates[79] += fwd_rates[1311];

  //rxn 1312
  //sp 59
  sp_rates[58] += fwd_rates[1312];
  //sp 60
  sp_rates[59] -= fwd_rates[1312];
  //sp 78
  sp_rates[77] -= fwd_rates[1312];
  //sp 80
  sp_rates[79] += fwd_rates[1312];

  //rxn 1313
  //sp 62
  sp_rates[61] += fwd_rates[1313];
  //sp 78
  sp_rates[77] -= fwd_rates[1313];
  //sp 80
  sp_rates[79] += fwd_rates[1313];
  //sp 64
  sp_rates[63] -= fwd_rates[1313];

  //rxn 1314
  //sp 65
  sp_rates[64] -= fwd_rates[1314];
  //sp 62
  sp_rates[61] += fwd_rates[1314];
  //sp 78
  sp_rates[77] -= fwd_rates[1314];
  //sp 80
  sp_rates[79] += fwd_rates[1314];

  //rxn 1315
  //sp 102
  sp_rates[101] += fwd_rates[1315];
  //sp 78
  sp_rates[77] -= fwd_rates[1315];
  //sp 103
  sp_rates[102] -= fwd_rates[1315];
  //sp 80
  sp_rates[79] += fwd_rates[1315];

  //rxn 1316
  //sp 97
  sp_rates[96] += fwd_rates[1316];
  //sp 100
  sp_rates[99] -= fwd_rates[1316];
  //sp 78
  sp_rates[77] -= fwd_rates[1316];
  //sp 80
  sp_rates[79] += fwd_rates[1316];

  //rxn 1317
  //sp 98
  sp_rates[97] += fwd_rates[1317];
  //sp 99
  sp_rates[98] -= fwd_rates[1317];
  //sp 78
  sp_rates[77] -= fwd_rates[1317];
  //sp 80
  sp_rates[79] += fwd_rates[1317];

  //rxn 1318
  //sp 94
  sp_rates[93] += fwd_rates[1318];
  //sp 78
  sp_rates[77] -= fwd_rates[1318];
  //sp 95
  sp_rates[94] -= fwd_rates[1318];
  //sp 80
  sp_rates[79] += fwd_rates[1318];

  //rxn 1319
  //sp 105
  sp_rates[104] -= fwd_rates[1319];
  //sp 104
  sp_rates[103] += fwd_rates[1319];
  //sp 78
  sp_rates[77] -= fwd_rates[1319];
  //sp 80
  sp_rates[79] += fwd_rates[1319];

  //rxn 1320
  //sp 107
  sp_rates[106] += fwd_rates[1320];
  //sp 109
  sp_rates[108] -= fwd_rates[1320];
  //sp 78
  sp_rates[77] -= fwd_rates[1320];
  //sp 80
  sp_rates[79] += fwd_rates[1320];

  //rxn 1321
  //sp 110
  sp_rates[109] += fwd_rates[1321];
  //sp 78
  sp_rates[77] -= fwd_rates[1321];
  //sp 111
  sp_rates[110] -= fwd_rates[1321];
  //sp 80
  sp_rates[79] += fwd_rates[1321];

  //rxn 1322
  //sp 22
  sp_rates[21] += fwd_rates[1322];
  //sp 78
  sp_rates[77] -= fwd_rates[1322];
  //sp 23
  sp_rates[22] -= fwd_rates[1322];
  //sp 80
  sp_rates[79] += fwd_rates[1322];

  //rxn 1323
  //sp 22
  sp_rates[21] += fwd_rates[1323];
  //sp 78
  sp_rates[77] -= fwd_rates[1323];
  //sp 80
  sp_rates[79] += fwd_rates[1323];
  //sp 24
  sp_rates[23] -= fwd_rates[1323];

  //rxn 1324
  //sp 25
  sp_rates[24] += fwd_rates[1324];
  //sp 26
  sp_rates[25] -= fwd_rates[1324];
  //sp 78
  sp_rates[77] -= fwd_rates[1324];
  //sp 80
  sp_rates[79] += fwd_rates[1324];

  //rxn 1325
  //sp 78
  sp_rates[77] -= fwd_rates[1325];
  //sp 46
  sp_rates[45] -= fwd_rates[1325];
  //sp 45
  sp_rates[44] += fwd_rates[1325];
  //sp 80
  sp_rates[79] += fwd_rates[1325];

  //rxn 1326
  //sp 45
  sp_rates[44] += fwd_rates[1326];
  //sp 78
  sp_rates[77] -= fwd_rates[1326];
  //sp 47
  sp_rates[46] -= fwd_rates[1326];
  //sp 80
  sp_rates[79] += fwd_rates[1326];

  //rxn 1327
  //sp 49
  sp_rates[48] -= fwd_rates[1327];
  //sp 78
  sp_rates[77] -= fwd_rates[1327];
  //sp 80
  sp_rates[79] += fwd_rates[1327];
  //sp 48
  sp_rates[47] += fwd_rates[1327];

  //rxn 1328
  //sp 111
  sp_rates[110] -= fwd_rates[1328];
  //sp 110
  sp_rates[109] += fwd_rates[1328];
  //sp 7
  sp_rates[6] -= fwd_rates[1328];
  //sp 8
  sp_rates[7] += fwd_rates[1328];

  //rxn 1329
  //sp 49
  sp_rates[48] -= fwd_rates[1329];
  //sp 48
  sp_rates[47] += fwd_rates[1329];
  //sp 7
  sp_rates[6] -= fwd_rates[1329];
  //sp 8
  sp_rates[7] += fwd_rates[1329];

  //rxn 1330
  //sp 9
  sp_rates[8] -= fwd_rates[1330];
  //sp 10
  sp_rates[9] += fwd_rates[1330];
  //sp 35
  sp_rates[34] -= fwd_rates[1330];
  //sp 34
  sp_rates[33] += fwd_rates[1330];

  //rxn 1331
  //sp 9
  sp_rates[8] -= fwd_rates[1331];
  //sp 10
  sp_rates[9] += fwd_rates[1331];
  //sp 105
  sp_rates[104] -= fwd_rates[1331];
  //sp 104
  sp_rates[103] += fwd_rates[1331];

  //rxn 1332
  //sp 9
  sp_rates[8] -= fwd_rates[1332];
  //sp 10
  sp_rates[9] += fwd_rates[1332];
  //sp 110
  sp_rates[109] += fwd_rates[1332];
  //sp 111
  sp_rates[110] -= fwd_rates[1332];

  //rxn 1333
  //sp 9
  sp_rates[8] -= fwd_rates[1333];
  //sp 10
  sp_rates[9] += fwd_rates[1333];
  //sp 45
  sp_rates[44] += fwd_rates[1333];
  //sp 47
  sp_rates[46] -= fwd_rates[1333];

  //rxn 1334
  //sp 9
  sp_rates[8] -= fwd_rates[1334];
  //sp 10
  sp_rates[9] += fwd_rates[1334];
  //sp 49
  sp_rates[48] -= fwd_rates[1334];
  //sp 48
  sp_rates[47] += fwd_rates[1334];

  //rxn 1335
  //sp 25
  sp_rates[24] -= fwd_rates[1335];
  //sp 26
  sp_rates[25] += fwd_rates[1335];
  //sp 30
  sp_rates[29] += fwd_rates[1335];
  //sp 31
  sp_rates[30] -= fwd_rates[1335];

  //rxn 1336
  //sp 25
  sp_rates[24] -= fwd_rates[1336];
  //sp 26
  sp_rates[25] += fwd_rates[1336];
  //sp 59
  sp_rates[58] += fwd_rates[1336];
  //sp 61
  sp_rates[60] -= fwd_rates[1336];

  //rxn 1337
  //sp 25
  sp_rates[24] -= fwd_rates[1337];
  //sp 26
  sp_rates[25] += fwd_rates[1337];
  //sp 59
  sp_rates[58] += fwd_rates[1337];
  //sp 60
  sp_rates[59] -= fwd_rates[1337];

  //rxn 1338
  //sp 25
  sp_rates[24] -= fwd_rates[1338];
  //sp 26
  sp_rates[25] += fwd_rates[1338];
  //sp 62
  sp_rates[61] += fwd_rates[1338];
  //sp 63
  sp_rates[62] -= fwd_rates[1338];

  //rxn 1339
  //sp 25
  sp_rates[24] -= fwd_rates[1339];
  //sp 26
  sp_rates[25] += fwd_rates[1339];
  //sp 102
  sp_rates[101] += fwd_rates[1339];
  //sp 103
  sp_rates[102] -= fwd_rates[1339];

  //rxn 1340
  //sp 25
  sp_rates[24] -= fwd_rates[1340];
  //sp 26
  sp_rates[25] += fwd_rates[1340];
  //sp 100
  sp_rates[99] -= fwd_rates[1340];
  //sp 97
  sp_rates[96] += fwd_rates[1340];

  //rxn 1341
  //sp 25
  sp_rates[24] -= fwd_rates[1341];
  //sp 26
  sp_rates[25] += fwd_rates[1341];
  //sp 99
  sp_rates[98] -= fwd_rates[1341];
  //sp 98
  sp_rates[97] += fwd_rates[1341];

  //rxn 1342
  //sp 25
  sp_rates[24] -= fwd_rates[1342];
  //sp 26
  sp_rates[25] += fwd_rates[1342];
  //sp 94
  sp_rates[93] += fwd_rates[1342];
  //sp 95
  sp_rates[94] -= fwd_rates[1342];

  //rxn 1343
  //sp 25
  sp_rates[24] -= fwd_rates[1343];
  //sp 26
  sp_rates[25] += fwd_rates[1343];
  //sp 105
  sp_rates[104] -= fwd_rates[1343];
  //sp 104
  sp_rates[103] += fwd_rates[1343];

  //rxn 1344
  //sp 25
  sp_rates[24] -= fwd_rates[1344];
  //sp 26
  sp_rates[25] += fwd_rates[1344];
  //sp 110
  sp_rates[109] += fwd_rates[1344];
  //sp 111
  sp_rates[110] -= fwd_rates[1344];

  //rxn 1345
  //sp 25
  sp_rates[24] -= fwd_rates[1345];
  //sp 26
  sp_rates[25] += fwd_rates[1345];
  //sp 45
  sp_rates[44] += fwd_rates[1345];
  //sp 46
  sp_rates[45] -= fwd_rates[1345];

  //rxn 1346
  //sp 25
  sp_rates[24] -= fwd_rates[1346];
  //sp 26
  sp_rates[25] += fwd_rates[1346];
  //sp 45
  sp_rates[44] += fwd_rates[1346];
  //sp 47
  sp_rates[46] -= fwd_rates[1346];

  //rxn 1347
  //sp 25
  sp_rates[24] -= fwd_rates[1347];
  //sp 26
  sp_rates[25] += fwd_rates[1347];
  //sp 79
  sp_rates[78] += fwd_rates[1347];
  //sp 80
  sp_rates[79] -= fwd_rates[1347];

  //rxn 1348
  //sp 25
  sp_rates[24] -= fwd_rates[1348];
  //sp 26
  sp_rates[25] += fwd_rates[1348];
  //sp 49
  sp_rates[48] -= fwd_rates[1348];
  //sp 48
  sp_rates[47] += fwd_rates[1348];

  //rxn 1349
  //sp 62
  sp_rates[61] += fwd_rates[1349];
  //sp 22
  sp_rates[21] -= fwd_rates[1349];
  //sp 23
  sp_rates[22] += fwd_rates[1349];
  //sp 64
  sp_rates[63] -= fwd_rates[1349];

  //rxn 1350
  //sp 65
  sp_rates[64] -= fwd_rates[1350];
  //sp 62
  sp_rates[61] += fwd_rates[1350];
  //sp 22
  sp_rates[21] -= fwd_rates[1350];
  //sp 23
  sp_rates[22] += fwd_rates[1350];

  //rxn 1351
  //sp 102
  sp_rates[101] += fwd_rates[1351];
  //sp 22
  sp_rates[21] -= fwd_rates[1351];
  //sp 103
  sp_rates[102] -= fwd_rates[1351];
  //sp 23
  sp_rates[22] += fwd_rates[1351];

  //rxn 1352
  //sp 105
  sp_rates[104] -= fwd_rates[1352];
  //sp 22
  sp_rates[21] -= fwd_rates[1352];
  //sp 23
  sp_rates[22] += fwd_rates[1352];
  //sp 104
  sp_rates[103] += fwd_rates[1352];

  //rxn 1353
  //sp 107
  sp_rates[106] += fwd_rates[1353];
  //sp 109
  sp_rates[108] -= fwd_rates[1353];
  //sp 22
  sp_rates[21] -= fwd_rates[1353];
  //sp 23
  sp_rates[22] += fwd_rates[1353];

  //rxn 1354
  //sp 110
  sp_rates[109] += fwd_rates[1354];
  //sp 22
  sp_rates[21] -= fwd_rates[1354];
  //sp 111
  sp_rates[110] -= fwd_rates[1354];
  //sp 23
  sp_rates[22] += fwd_rates[1354];

  //rxn 1355
  //sp 46
  sp_rates[45] -= fwd_rates[1355];
  //sp 45
  sp_rates[44] += fwd_rates[1355];
  //sp 22
  sp_rates[21] -= fwd_rates[1355];
  //sp 23
  sp_rates[22] += fwd_rates[1355];

  //rxn 1356
  //sp 45
  sp_rates[44] += fwd_rates[1356];
  //sp 22
  sp_rates[21] -= fwd_rates[1356];
  //sp 47
  sp_rates[46] -= fwd_rates[1356];
  //sp 23
  sp_rates[22] += fwd_rates[1356];

  //rxn 1357
  //sp 79
  sp_rates[78] += fwd_rates[1357];
  //sp 22
  sp_rates[21] -= fwd_rates[1357];
  //sp 23
  sp_rates[22] += fwd_rates[1357];
  //sp 80
  sp_rates[79] -= fwd_rates[1357];

  //rxn 1358
  //sp 49
  sp_rates[48] -= fwd_rates[1358];
  //sp 22
  sp_rates[21] -= fwd_rates[1358];
  //sp 23
  sp_rates[22] += fwd_rates[1358];
  //sp 48
  sp_rates[47] += fwd_rates[1358];

  //rxn 1359
  //sp 34
  sp_rates[33] += fwd_rates[1359];
  //sp 35
  sp_rates[34] -= fwd_rates[1359];
  //sp 22
  sp_rates[21] -= fwd_rates[1359];
  //sp 24
  sp_rates[23] += fwd_rates[1359];

  //rxn 1360
  //sp 30
  sp_rates[29] += fwd_rates[1360];
  //sp 22
  sp_rates[21] -= fwd_rates[1360];
  //sp 31
  sp_rates[30] -= fwd_rates[1360];
  //sp 24
  sp_rates[23] += fwd_rates[1360];

  //rxn 1361
  //sp 59
  sp_rates[58] += fwd_rates[1361];
  //sp 61
  sp_rates[60] -= fwd_rates[1361];
  //sp 22
  sp_rates[21] -= fwd_rates[1361];
  //sp 24
  sp_rates[23] += fwd_rates[1361];

  //rxn 1362
  //sp 59
  sp_rates[58] += fwd_rates[1362];
  //sp 60
  sp_rates[59] -= fwd_rates[1362];
  //sp 22
  sp_rates[21] -= fwd_rates[1362];
  //sp 24
  sp_rates[23] += fwd_rates[1362];

  //rxn 1363
  //sp 62
  sp_rates[61] += fwd_rates[1363];
  //sp 24
  sp_rates[23] += fwd_rates[1363];
  //sp 22
  sp_rates[21] -= fwd_rates[1363];
  //sp 64
  sp_rates[63] -= fwd_rates[1363];

  //rxn 1364
  //sp 65
  sp_rates[64] -= fwd_rates[1364];
  //sp 62
  sp_rates[61] += fwd_rates[1364];
  //sp 22
  sp_rates[21] -= fwd_rates[1364];
  //sp 24
  sp_rates[23] += fwd_rates[1364];

  //rxn 1365
  //sp 62
  sp_rates[61] += fwd_rates[1365];
  //sp 22
  sp_rates[21] -= fwd_rates[1365];
  //sp 63
  sp_rates[62] -= fwd_rates[1365];
  //sp 24
  sp_rates[23] += fwd_rates[1365];

  //rxn 1366
  //sp 102
  sp_rates[101] += fwd_rates[1366];
  //sp 22
  sp_rates[21] -= fwd_rates[1366];
  //sp 103
  sp_rates[102] -= fwd_rates[1366];
  //sp 24
  sp_rates[23] += fwd_rates[1366];

  //rxn 1367
  //sp 97
  sp_rates[96] += fwd_rates[1367];
  //sp 100
  sp_rates[99] -= fwd_rates[1367];
  //sp 22
  sp_rates[21] -= fwd_rates[1367];
  //sp 24
  sp_rates[23] += fwd_rates[1367];

  //rxn 1368
  //sp 98
  sp_rates[97] += fwd_rates[1368];
  //sp 99
  sp_rates[98] -= fwd_rates[1368];
  //sp 22
  sp_rates[21] -= fwd_rates[1368];
  //sp 24
  sp_rates[23] += fwd_rates[1368];

  //rxn 1369
  //sp 94
  sp_rates[93] += fwd_rates[1369];
  //sp 22
  sp_rates[21] -= fwd_rates[1369];
  //sp 95
  sp_rates[94] -= fwd_rates[1369];
  //sp 24
  sp_rates[23] += fwd_rates[1369];

  //rxn 1370
  //sp 105
  sp_rates[104] -= fwd_rates[1370];
  //sp 104
  sp_rates[103] += fwd_rates[1370];
  //sp 22
  sp_rates[21] -= fwd_rates[1370];
  //sp 24
  sp_rates[23] += fwd_rates[1370];

  //rxn 1371
  //sp 107
  sp_rates[106] += fwd_rates[1371];
  //sp 109
  sp_rates[108] -= fwd_rates[1371];
  //sp 22
  sp_rates[21] -= fwd_rates[1371];
  //sp 24
  sp_rates[23] += fwd_rates[1371];

  //rxn 1372
  //sp 110
  sp_rates[109] += fwd_rates[1372];
  //sp 22
  sp_rates[21] -= fwd_rates[1372];
  //sp 111
  sp_rates[110] -= fwd_rates[1372];
  //sp 24
  sp_rates[23] += fwd_rates[1372];

  //rxn 1373
  //sp 46
  sp_rates[45] -= fwd_rates[1373];
  //sp 45
  sp_rates[44] += fwd_rates[1373];
  //sp 22
  sp_rates[21] -= fwd_rates[1373];
  //sp 24
  sp_rates[23] += fwd_rates[1373];

  //rxn 1374
  //sp 45
  sp_rates[44] += fwd_rates[1374];
  //sp 22
  sp_rates[21] -= fwd_rates[1374];
  //sp 47
  sp_rates[46] -= fwd_rates[1374];
  //sp 24
  sp_rates[23] += fwd_rates[1374];

  //rxn 1375
  //sp 24
  sp_rates[23] += fwd_rates[1375];
  //sp 22
  sp_rates[21] -= fwd_rates[1375];
  //sp 79
  sp_rates[78] += fwd_rates[1375];
  //sp 80
  sp_rates[79] -= fwd_rates[1375];

  //rxn 1376
  //sp 49
  sp_rates[48] -= fwd_rates[1376];
  //sp 48
  sp_rates[47] += fwd_rates[1376];
  //sp 22
  sp_rates[21] -= fwd_rates[1376];
  //sp 24
  sp_rates[23] += fwd_rates[1376];

  //rxn 1377
  //sp 46
  sp_rates[45] += fwd_rates[1377];
  //sp 45
  sp_rates[44] -= fwd_rates[1377];
  //sp 14
  sp_rates[13] += fwd_rates[1377];
  //sp 15
  sp_rates[14] -= fwd_rates[1377];

  //rxn 1378
  //sp 34
  sp_rates[33] += fwd_rates[1378];
  //sp 35
  sp_rates[34] -= fwd_rates[1378];
  //sp 45
  sp_rates[44] -= fwd_rates[1378];
  //sp 46
  sp_rates[45] += fwd_rates[1378];

  //rxn 1379
  //sp 46
  sp_rates[45] += fwd_rates[1379];
  //sp 45
  sp_rates[44] -= fwd_rates[1379];
  //sp 30
  sp_rates[29] += fwd_rates[1379];
  //sp 31
  sp_rates[30] -= fwd_rates[1379];

  //rxn 1380
  //sp 59
  sp_rates[58] += fwd_rates[1380];
  //sp 45
  sp_rates[44] -= fwd_rates[1380];
  //sp 46
  sp_rates[45] += fwd_rates[1380];
  //sp 61
  sp_rates[60] -= fwd_rates[1380];

  //rxn 1381
  //sp 59
  sp_rates[58] += fwd_rates[1381];
  //sp 60
  sp_rates[59] -= fwd_rates[1381];
  //sp 45
  sp_rates[44] -= fwd_rates[1381];
  //sp 46
  sp_rates[45] += fwd_rates[1381];

  //rxn 1382
  //sp 62
  sp_rates[61] += fwd_rates[1382];
  //sp 45
  sp_rates[44] -= fwd_rates[1382];
  //sp 46
  sp_rates[45] += fwd_rates[1382];
  //sp 64
  sp_rates[63] -= fwd_rates[1382];

  //rxn 1383
  //sp 65
  sp_rates[64] -= fwd_rates[1383];
  //sp 62
  sp_rates[61] += fwd_rates[1383];
  //sp 45
  sp_rates[44] -= fwd_rates[1383];
  //sp 46
  sp_rates[45] += fwd_rates[1383];

  //rxn 1384
  //sp 62
  sp_rates[61] += fwd_rates[1384];
  //sp 45
  sp_rates[44] -= fwd_rates[1384];
  //sp 46
  sp_rates[45] += fwd_rates[1384];
  //sp 63
  sp_rates[62] -= fwd_rates[1384];

  //rxn 1385
  //sp 102
  sp_rates[101] += fwd_rates[1385];
  //sp 45
  sp_rates[44] -= fwd_rates[1385];
  //sp 46
  sp_rates[45] += fwd_rates[1385];
  //sp 103
  sp_rates[102] -= fwd_rates[1385];

  //rxn 1386
  //sp 97
  sp_rates[96] += fwd_rates[1386];
  //sp 100
  sp_rates[99] -= fwd_rates[1386];
  //sp 45
  sp_rates[44] -= fwd_rates[1386];
  //sp 46
  sp_rates[45] += fwd_rates[1386];

  //rxn 1387
  //sp 98
  sp_rates[97] += fwd_rates[1387];
  //sp 99
  sp_rates[98] -= fwd_rates[1387];
  //sp 45
  sp_rates[44] -= fwd_rates[1387];
  //sp 46
  sp_rates[45] += fwd_rates[1387];

  //rxn 1388
  //sp 94
  sp_rates[93] += fwd_rates[1388];
  //sp 45
  sp_rates[44] -= fwd_rates[1388];
  //sp 46
  sp_rates[45] += fwd_rates[1388];
  //sp 95
  sp_rates[94] -= fwd_rates[1388];

  //rxn 1389
  //sp 105
  sp_rates[104] -= fwd_rates[1389];
  //sp 45
  sp_rates[44] -= fwd_rates[1389];
  //sp 46
  sp_rates[45] += fwd_rates[1389];
  //sp 104
  sp_rates[103] += fwd_rates[1389];

  //rxn 1390
  //sp 109
  sp_rates[108] -= fwd_rates[1390];
  //sp 107
  sp_rates[106] += fwd_rates[1390];
  //sp 45
  sp_rates[44] -= fwd_rates[1390];
  //sp 46
  sp_rates[45] += fwd_rates[1390];

  //rxn 1391
  //sp 110
  sp_rates[109] += fwd_rates[1391];
  //sp 45
  sp_rates[44] -= fwd_rates[1391];
  //sp 46
  sp_rates[45] += fwd_rates[1391];
  //sp 111
  sp_rates[110] -= fwd_rates[1391];

  //rxn 1392
  //sp 46
  sp_rates[45] += fwd_rates[1392];
  //sp 45
  sp_rates[44] -= fwd_rates[1392];
  //sp 22
  sp_rates[21] += fwd_rates[1392];
  //sp 23
  sp_rates[22] -= fwd_rates[1392];

  //rxn 1393
  //sp 46
  sp_rates[45] += fwd_rates[1393];
  //sp 45
  sp_rates[44] -= fwd_rates[1393];
  //sp 22
  sp_rates[21] += fwd_rates[1393];
  //sp 24
  sp_rates[23] -= fwd_rates[1393];

  //rxn 1394
  //sp 25
  sp_rates[24] += fwd_rates[1394];
  //sp 26
  sp_rates[25] -= fwd_rates[1394];
  //sp 45
  sp_rates[44] -= fwd_rates[1394];
  //sp 46
  sp_rates[45] += fwd_rates[1394];

  //rxn 1395
  //sp 46
  sp_rates[45] += fwd_rates[1395];
  //sp 47
  sp_rates[46] -= fwd_rates[1395];

  //rxn 1396
  //sp 45
  sp_rates[44] -= fwd_rates[1396];
  //sp 46
  sp_rates[45] += fwd_rates[1396];
  //sp 79
  sp_rates[78] += fwd_rates[1396];
  //sp 80
  sp_rates[79] -= fwd_rates[1396];

  //rxn 1397
  //sp 49
  sp_rates[48] -= fwd_rates[1397];
  //sp 45
  sp_rates[44] -= fwd_rates[1397];
  //sp 46
  sp_rates[45] += fwd_rates[1397];
  //sp 48
  sp_rates[47] += fwd_rates[1397];

  //rxn 1398
  //sp 47
  sp_rates[46] += fwd_rates[1398];
  //sp 45
  sp_rates[44] -= fwd_rates[1398];
  //sp 14
  sp_rates[13] += fwd_rates[1398];
  //sp 15
  sp_rates[14] -= fwd_rates[1398];

  //rxn 1399
  //sp 34
  sp_rates[33] += fwd_rates[1399];
  //sp 35
  sp_rates[34] -= fwd_rates[1399];
  //sp 45
  sp_rates[44] -= fwd_rates[1399];
  //sp 47
  sp_rates[46] += fwd_rates[1399];

  //rxn 1400
  //sp 47
  sp_rates[46] += fwd_rates[1400];
  //sp 45
  sp_rates[44] -= fwd_rates[1400];
  //sp 30
  sp_rates[29] += fwd_rates[1400];
  //sp 31
  sp_rates[30] -= fwd_rates[1400];

  //rxn 1401
  //sp 47
  sp_rates[46] += fwd_rates[1401];
  //sp 59
  sp_rates[58] += fwd_rates[1401];
  //sp 45
  sp_rates[44] -= fwd_rates[1401];
  //sp 61
  sp_rates[60] -= fwd_rates[1401];

  //rxn 1402
  //sp 59
  sp_rates[58] += fwd_rates[1402];
  //sp 60
  sp_rates[59] -= fwd_rates[1402];
  //sp 45
  sp_rates[44] -= fwd_rates[1402];
  //sp 47
  sp_rates[46] += fwd_rates[1402];

  //rxn 1403
  //sp 45
  sp_rates[44] -= fwd_rates[1403];
  //sp 62
  sp_rates[61] += fwd_rates[1403];
  //sp 47
  sp_rates[46] += fwd_rates[1403];
  //sp 64
  sp_rates[63] -= fwd_rates[1403];

  //rxn 1404
  //sp 65
  sp_rates[64] -= fwd_rates[1404];
  //sp 45
  sp_rates[44] -= fwd_rates[1404];
  //sp 62
  sp_rates[61] += fwd_rates[1404];
  //sp 47
  sp_rates[46] += fwd_rates[1404];

  //rxn 1405
  //sp 47
  sp_rates[46] += fwd_rates[1405];
  //sp 45
  sp_rates[44] -= fwd_rates[1405];
  //sp 62
  sp_rates[61] += fwd_rates[1405];
  //sp 63
  sp_rates[62] -= fwd_rates[1405];

  //rxn 1406
  //sp 47
  sp_rates[46] += fwd_rates[1406];
  //sp 45
  sp_rates[44] -= fwd_rates[1406];
  //sp 102
  sp_rates[101] += fwd_rates[1406];
  //sp 103
  sp_rates[102] -= fwd_rates[1406];

  //rxn 1407
  //sp 97
  sp_rates[96] += fwd_rates[1407];
  //sp 100
  sp_rates[99] -= fwd_rates[1407];
  //sp 45
  sp_rates[44] -= fwd_rates[1407];
  //sp 47
  sp_rates[46] += fwd_rates[1407];

  //rxn 1408
  //sp 98
  sp_rates[97] += fwd_rates[1408];
  //sp 99
  sp_rates[98] -= fwd_rates[1408];
  //sp 45
  sp_rates[44] -= fwd_rates[1408];
  //sp 47
  sp_rates[46] += fwd_rates[1408];

  //rxn 1409
  //sp 47
  sp_rates[46] += fwd_rates[1409];
  //sp 45
  sp_rates[44] -= fwd_rates[1409];
  //sp 94
  sp_rates[93] += fwd_rates[1409];
  //sp 95
  sp_rates[94] -= fwd_rates[1409];

  //rxn 1410
  //sp 105
  sp_rates[104] -= fwd_rates[1410];
  //sp 45
  sp_rates[44] -= fwd_rates[1410];
  //sp 47
  sp_rates[46] += fwd_rates[1410];
  //sp 104
  sp_rates[103] += fwd_rates[1410];

  //rxn 1411
  //sp 109
  sp_rates[108] -= fwd_rates[1411];
  //sp 107
  sp_rates[106] += fwd_rates[1411];
  //sp 45
  sp_rates[44] -= fwd_rates[1411];
  //sp 47
  sp_rates[46] += fwd_rates[1411];

  //rxn 1412
  //sp 47
  sp_rates[46] += fwd_rates[1412];
  //sp 45
  sp_rates[44] -= fwd_rates[1412];
  //sp 110
  sp_rates[109] += fwd_rates[1412];
  //sp 111
  sp_rates[110] -= fwd_rates[1412];

  //rxn 1413
  //sp 45
  sp_rates[44] -= fwd_rates[1413];
  //sp 6
  sp_rates[5] -= fwd_rates[1413];
  //sp 47
  sp_rates[46] += fwd_rates[1413];
  //sp 8
  sp_rates[7] += fwd_rates[1413];

  //rxn 1414
  //sp 9
  sp_rates[8] += fwd_rates[1414];
  //sp 10
  sp_rates[9] -= fwd_rates[1414];
  //sp 45
  sp_rates[44] -= fwd_rates[1414];
  //sp 47
  sp_rates[46] += fwd_rates[1414];

  //rxn 1415
  //sp 47
  sp_rates[46] += fwd_rates[1415];
  //sp 45
  sp_rates[44] -= fwd_rates[1415];
  //sp 22
  sp_rates[21] += fwd_rates[1415];
  //sp 23
  sp_rates[22] -= fwd_rates[1415];

  //rxn 1416
  //sp 45
  sp_rates[44] -= fwd_rates[1416];
  //sp 22
  sp_rates[21] += fwd_rates[1416];
  //sp 47
  sp_rates[46] += fwd_rates[1416];
  //sp 24
  sp_rates[23] -= fwd_rates[1416];

  //rxn 1417
  //sp 25
  sp_rates[24] += fwd_rates[1417];
  //sp 26
  sp_rates[25] -= fwd_rates[1417];
  //sp 45
  sp_rates[44] -= fwd_rates[1417];
  //sp 47
  sp_rates[46] += fwd_rates[1417];

  //rxn 1418
  //sp 46
  sp_rates[45] -= fwd_rates[1418];
  //sp 47
  sp_rates[46] += fwd_rates[1418];

  //rxn 1419
  //sp 79
  sp_rates[78] += fwd_rates[1419];
  //sp 45
  sp_rates[44] -= fwd_rates[1419];
  //sp 47
  sp_rates[46] += fwd_rates[1419];
  //sp 80
  sp_rates[79] -= fwd_rates[1419];

  //rxn 1420
  //sp 49
  sp_rates[48] -= fwd_rates[1420];
  //sp 45
  sp_rates[44] -= fwd_rates[1420];
  //sp 47
  sp_rates[46] += fwd_rates[1420];
  //sp 48
  sp_rates[47] += fwd_rates[1420];

  //rxn 1421
  //sp 84
  sp_rates[83] += fwd_rates[1421];
  //sp 82
  sp_rates[81] -= fwd_rates[1421];
  //sp 3
  sp_rates[2] += fwd_rates[1421];
  //sp 4
  sp_rates[3] -= fwd_rates[1421];

  //rxn 1422
  //sp 82
  sp_rates[81] -= fwd_rates[1422];
  //sp 84
  sp_rates[83] += fwd_rates[1422];
  //sp 14
  sp_rates[13] += fwd_rates[1422];
  //sp 15
  sp_rates[14] -= fwd_rates[1422];

  //rxn 1423
  //sp 82
  sp_rates[81] -= fwd_rates[1423];
  //sp 35
  sp_rates[34] -= fwd_rates[1423];
  //sp 84
  sp_rates[83] += fwd_rates[1423];
  //sp 34
  sp_rates[33] += fwd_rates[1423];

  //rxn 1424
  //sp 82
  sp_rates[81] -= fwd_rates[1424];
  //sp 84
  sp_rates[83] += fwd_rates[1424];
  //sp 30
  sp_rates[29] += fwd_rates[1424];
  //sp 31
  sp_rates[30] -= fwd_rates[1424];

  //rxn 1425
  //sp 82
  sp_rates[81] -= fwd_rates[1425];
  //sp 59
  sp_rates[58] += fwd_rates[1425];
  //sp 84
  sp_rates[83] += fwd_rates[1425];
  //sp 61
  sp_rates[60] -= fwd_rates[1425];

  //rxn 1426
  //sp 84
  sp_rates[83] += fwd_rates[1426];
  //sp 82
  sp_rates[81] -= fwd_rates[1426];
  //sp 59
  sp_rates[58] += fwd_rates[1426];
  //sp 60
  sp_rates[59] -= fwd_rates[1426];

  //rxn 1427
  //sp 82
  sp_rates[81] -= fwd_rates[1427];
  //sp 84
  sp_rates[83] += fwd_rates[1427];
  //sp 62
  sp_rates[61] += fwd_rates[1427];
  //sp 64
  sp_rates[63] -= fwd_rates[1427];

  //rxn 1428
  //sp 65
  sp_rates[64] -= fwd_rates[1428];
  //sp 82
  sp_rates[81] -= fwd_rates[1428];
  //sp 84
  sp_rates[83] += fwd_rates[1428];
  //sp 62
  sp_rates[61] += fwd_rates[1428];

  //rxn 1429
  //sp 82
  sp_rates[81] -= fwd_rates[1429];
  //sp 84
  sp_rates[83] += fwd_rates[1429];
  //sp 62
  sp_rates[61] += fwd_rates[1429];
  //sp 63
  sp_rates[62] -= fwd_rates[1429];

  //rxn 1430
  //sp 82
  sp_rates[81] -= fwd_rates[1430];
  //sp 84
  sp_rates[83] += fwd_rates[1430];
  //sp 102
  sp_rates[101] += fwd_rates[1430];
  //sp 103
  sp_rates[102] -= fwd_rates[1430];

  //rxn 1431
  //sp 97
  sp_rates[96] += fwd_rates[1431];
  //sp 82
  sp_rates[81] -= fwd_rates[1431];
  //sp 84
  sp_rates[83] += fwd_rates[1431];
  //sp 100
  sp_rates[99] -= fwd_rates[1431];

  //rxn 1432
  //sp 82
  sp_rates[81] -= fwd_rates[1432];
  //sp 99
  sp_rates[98] -= fwd_rates[1432];
  //sp 84
  sp_rates[83] += fwd_rates[1432];
  //sp 98
  sp_rates[97] += fwd_rates[1432];

  //rxn 1433
  //sp 82
  sp_rates[81] -= fwd_rates[1433];
  //sp 84
  sp_rates[83] += fwd_rates[1433];
  //sp 94
  sp_rates[93] += fwd_rates[1433];
  //sp 95
  sp_rates[94] -= fwd_rates[1433];

  //rxn 1434
  //sp 105
  sp_rates[104] -= fwd_rates[1434];
  //sp 82
  sp_rates[81] -= fwd_rates[1434];
  //sp 84
  sp_rates[83] += fwd_rates[1434];
  //sp 104
  sp_rates[103] += fwd_rates[1434];

  //rxn 1435
  //sp 82
  sp_rates[81] -= fwd_rates[1435];
  //sp 107
  sp_rates[106] += fwd_rates[1435];
  //sp 84
  sp_rates[83] += fwd_rates[1435];
  //sp 109
  sp_rates[108] -= fwd_rates[1435];

  //rxn 1436
  //sp 82
  sp_rates[81] -= fwd_rates[1436];
  //sp 84
  sp_rates[83] += fwd_rates[1436];
  //sp 110
  sp_rates[109] += fwd_rates[1436];
  //sp 111
  sp_rates[110] -= fwd_rates[1436];

  //rxn 1437
  //sp 82
  sp_rates[81] -= fwd_rates[1437];
  //sp 84
  sp_rates[83] += fwd_rates[1437];
  //sp 6
  sp_rates[5] -= fwd_rates[1437];
  //sp 8
  sp_rates[7] += fwd_rates[1437];

  //rxn 1438
  //sp 82
  sp_rates[81] -= fwd_rates[1438];
  //sp 84
  sp_rates[83] += fwd_rates[1438];
  //sp 5
  sp_rates[4] -= fwd_rates[1438];
  //sp 10
  sp_rates[9] += fwd_rates[1438];

  //rxn 1439
  //sp 82
  sp_rates[81] -= fwd_rates[1439];
  //sp 84
  sp_rates[83] += fwd_rates[1439];
  //sp 22
  sp_rates[21] += fwd_rates[1439];
  //sp 23
  sp_rates[22] -= fwd_rates[1439];

  //rxn 1440
  //sp 82
  sp_rates[81] -= fwd_rates[1440];
  //sp 84
  sp_rates[83] += fwd_rates[1440];
  //sp 22
  sp_rates[21] += fwd_rates[1440];
  //sp 24
  sp_rates[23] -= fwd_rates[1440];

  //rxn 1441
  //sp 82
  sp_rates[81] -= fwd_rates[1441];
  //sp 26
  sp_rates[25] -= fwd_rates[1441];
  //sp 84
  sp_rates[83] += fwd_rates[1441];
  //sp 25
  sp_rates[24] += fwd_rates[1441];

  //rxn 1442
  //sp 82
  sp_rates[81] -= fwd_rates[1442];
  //sp 84
  sp_rates[83] += fwd_rates[1442];
  //sp 45
  sp_rates[44] += fwd_rates[1442];
  //sp 46
  sp_rates[45] -= fwd_rates[1442];

  //rxn 1443
  //sp 82
  sp_rates[81] -= fwd_rates[1443];
  //sp 84
  sp_rates[83] += fwd_rates[1443];
  //sp 45
  sp_rates[44] += fwd_rates[1443];
  //sp 47
  sp_rates[46] -= fwd_rates[1443];

  //rxn 1444
  //sp 82
  sp_rates[81] -= fwd_rates[1444];
  //sp 84
  sp_rates[83] += fwd_rates[1444];
  //sp 79
  sp_rates[78] += fwd_rates[1444];
  //sp 80
  sp_rates[79] -= fwd_rates[1444];

  //rxn 1445
  //sp 49
  sp_rates[48] -= fwd_rates[1445];
  //sp 82
  sp_rates[81] -= fwd_rates[1445];
  //sp 84
  sp_rates[83] += fwd_rates[1445];
  //sp 48
  sp_rates[47] += fwd_rates[1445];

  //rxn 1446
  //sp 34
  sp_rates[33] += fwd_rates[1446];
  //sp 35
  sp_rates[34] -= fwd_rates[1446];
  //sp 11
  sp_rates[10] += fwd_rates[1446];
  //sp 82
  sp_rates[81] -= fwd_rates[1446];
  //sp 31
  sp_rates[30] += fwd_rates[1446];

  //rxn 1447
  //sp 82
  sp_rates[81] -= fwd_rates[1447];
  //sp 11
  sp_rates[10] += fwd_rates[1447];
  //sp 30
  sp_rates[29] += fwd_rates[1447];

  //rxn 1448
  //sp 11
  sp_rates[10] += fwd_rates[1448];
  //sp 82
  sp_rates[81] -= fwd_rates[1448];
  //sp 59
  sp_rates[58] += fwd_rates[1448];
  //sp 61
  sp_rates[60] -= fwd_rates[1448];
  //sp 31
  sp_rates[30] += fwd_rates[1448];

  //rxn 1449
  //sp 11
  sp_rates[10] += fwd_rates[1449];
  //sp 82
  sp_rates[81] -= fwd_rates[1449];
  //sp 59
  sp_rates[58] += fwd_rates[1449];
  //sp 60
  sp_rates[59] -= fwd_rates[1449];
  //sp 31
  sp_rates[30] += fwd_rates[1449];

  //rxn 1450
  //sp 11
  sp_rates[10] += fwd_rates[1450];
  //sp 82
  sp_rates[81] -= fwd_rates[1450];
  //sp 62
  sp_rates[61] += fwd_rates[1450];
  //sp 31
  sp_rates[30] += fwd_rates[1450];
  //sp 64
  sp_rates[63] -= fwd_rates[1450];

  //rxn 1451
  //sp 65
  sp_rates[64] -= fwd_rates[1451];
  //sp 11
  sp_rates[10] += fwd_rates[1451];
  //sp 82
  sp_rates[81] -= fwd_rates[1451];
  //sp 62
  sp_rates[61] += fwd_rates[1451];
  //sp 31
  sp_rates[30] += fwd_rates[1451];

  //rxn 1452
  //sp 11
  sp_rates[10] += fwd_rates[1452];
  //sp 82
  sp_rates[81] -= fwd_rates[1452];
  //sp 31
  sp_rates[30] += fwd_rates[1452];
  //sp 62
  sp_rates[61] += fwd_rates[1452];
  //sp 63
  sp_rates[62] -= fwd_rates[1452];

  //rxn 1453
  //sp 102
  sp_rates[101] += fwd_rates[1453];
  //sp 103
  sp_rates[102] -= fwd_rates[1453];
  //sp 11
  sp_rates[10] += fwd_rates[1453];
  //sp 82
  sp_rates[81] -= fwd_rates[1453];
  //sp 31
  sp_rates[30] += fwd_rates[1453];

  //rxn 1454
  //sp 97
  sp_rates[96] += fwd_rates[1454];
  //sp 100
  sp_rates[99] -= fwd_rates[1454];
  //sp 11
  sp_rates[10] += fwd_rates[1454];
  //sp 82
  sp_rates[81] -= fwd_rates[1454];
  //sp 31
  sp_rates[30] += fwd_rates[1454];

  //rxn 1455
  //sp 98
  sp_rates[97] += fwd_rates[1455];
  //sp 99
  sp_rates[98] -= fwd_rates[1455];
  //sp 11
  sp_rates[10] += fwd_rates[1455];
  //sp 82
  sp_rates[81] -= fwd_rates[1455];
  //sp 31
  sp_rates[30] += fwd_rates[1455];

  //rxn 1456
  //sp 11
  sp_rates[10] += fwd_rates[1456];
  //sp 82
  sp_rates[81] -= fwd_rates[1456];
  //sp 31
  sp_rates[30] += fwd_rates[1456];
  //sp 94
  sp_rates[93] += fwd_rates[1456];
  //sp 95
  sp_rates[94] -= fwd_rates[1456];

  //rxn 1457
  //sp 104
  sp_rates[103] += fwd_rates[1457];
  //sp 105
  sp_rates[104] -= fwd_rates[1457];
  //sp 11
  sp_rates[10] += fwd_rates[1457];
  //sp 82
  sp_rates[81] -= fwd_rates[1457];
  //sp 31
  sp_rates[30] += fwd_rates[1457];

  //rxn 1458
  //sp 11
  sp_rates[10] += fwd_rates[1458];
  //sp 107
  sp_rates[106] += fwd_rates[1458];
  //sp 109
  sp_rates[108] -= fwd_rates[1458];
  //sp 82
  sp_rates[81] -= fwd_rates[1458];
  //sp 31
  sp_rates[30] += fwd_rates[1458];

  //rxn 1459
  //sp 11
  sp_rates[10] += fwd_rates[1459];
  //sp 110
  sp_rates[109] += fwd_rates[1459];
  //sp 111
  sp_rates[110] -= fwd_rates[1459];
  //sp 82
  sp_rates[81] -= fwd_rates[1459];
  //sp 31
  sp_rates[30] += fwd_rates[1459];

  //rxn 1460
  //sp 6
  sp_rates[5] -= fwd_rates[1460];
  //sp 8
  sp_rates[7] += fwd_rates[1460];
  //sp 11
  sp_rates[10] += fwd_rates[1460];
  //sp 82
  sp_rates[81] -= fwd_rates[1460];
  //sp 31
  sp_rates[30] += fwd_rates[1460];

  //rxn 1461
  //sp 5
  sp_rates[4] -= fwd_rates[1461];
  //sp 10
  sp_rates[9] += fwd_rates[1461];
  //sp 11
  sp_rates[10] += fwd_rates[1461];
  //sp 82
  sp_rates[81] -= fwd_rates[1461];
  //sp 31
  sp_rates[30] += fwd_rates[1461];

  //rxn 1462
  //sp 11
  sp_rates[10] += fwd_rates[1462];
  //sp 82
  sp_rates[81] -= fwd_rates[1462];
  //sp 22
  sp_rates[21] += fwd_rates[1462];
  //sp 23
  sp_rates[22] -= fwd_rates[1462];
  //sp 31
  sp_rates[30] += fwd_rates[1462];

  //rxn 1463
  //sp 11
  sp_rates[10] += fwd_rates[1463];
  //sp 82
  sp_rates[81] -= fwd_rates[1463];
  //sp 22
  sp_rates[21] += fwd_rates[1463];
  //sp 24
  sp_rates[23] -= fwd_rates[1463];
  //sp 31
  sp_rates[30] += fwd_rates[1463];

  //rxn 1464
  //sp 11
  sp_rates[10] += fwd_rates[1464];
  //sp 82
  sp_rates[81] -= fwd_rates[1464];
  //sp 25
  sp_rates[24] += fwd_rates[1464];
  //sp 26
  sp_rates[25] -= fwd_rates[1464];
  //sp 31
  sp_rates[30] += fwd_rates[1464];

  //rxn 1465
  //sp 11
  sp_rates[10] += fwd_rates[1465];
  //sp 45
  sp_rates[44] += fwd_rates[1465];
  //sp 46
  sp_rates[45] -= fwd_rates[1465];
  //sp 82
  sp_rates[81] -= fwd_rates[1465];
  //sp 31
  sp_rates[30] += fwd_rates[1465];

  //rxn 1466
  //sp 11
  sp_rates[10] += fwd_rates[1466];
  //sp 45
  sp_rates[44] += fwd_rates[1466];
  //sp 47
  sp_rates[46] -= fwd_rates[1466];
  //sp 82
  sp_rates[81] -= fwd_rates[1466];
  //sp 31
  sp_rates[30] += fwd_rates[1466];

  //rxn 1467
  //sp 11
  sp_rates[10] += fwd_rates[1467];
  //sp 79
  sp_rates[78] += fwd_rates[1467];
  //sp 80
  sp_rates[79] -= fwd_rates[1467];
  //sp 82
  sp_rates[81] -= fwd_rates[1467];
  //sp 31
  sp_rates[30] += fwd_rates[1467];

  //rxn 1468
  //sp 11
  sp_rates[10] += fwd_rates[1468];
  //sp 48
  sp_rates[47] += fwd_rates[1468];
  //sp 49
  sp_rates[48] -= fwd_rates[1468];
  //sp 82
  sp_rates[81] -= fwd_rates[1468];
  //sp 31
  sp_rates[30] += fwd_rates[1468];

  //rxn 1469
  //sp 49
  sp_rates[48] += fwd_rates[1469];
  //sp 3
  sp_rates[2] += fwd_rates[1469];
  //sp 4
  sp_rates[3] -= fwd_rates[1469];
  //sp 48
  sp_rates[47] -= fwd_rates[1469];

  //rxn 1470
  //sp 49
  sp_rates[48] += fwd_rates[1470];
  //sp 14
  sp_rates[13] += fwd_rates[1470];
  //sp 15
  sp_rates[14] -= fwd_rates[1470];
  //sp 48
  sp_rates[47] -= fwd_rates[1470];

  //rxn 1471
  //sp 49
  sp_rates[48] += fwd_rates[1471];
  //sp 34
  sp_rates[33] += fwd_rates[1471];
  //sp 35
  sp_rates[34] -= fwd_rates[1471];
  //sp 48
  sp_rates[47] -= fwd_rates[1471];

  //rxn 1472
  //sp 49
  sp_rates[48] += fwd_rates[1472];
  //sp 30
  sp_rates[29] += fwd_rates[1472];
  //sp 31
  sp_rates[30] -= fwd_rates[1472];
  //sp 48
  sp_rates[47] -= fwd_rates[1472];

  //rxn 1473
  //sp 49
  sp_rates[48] += fwd_rates[1473];
  //sp 59
  sp_rates[58] += fwd_rates[1473];
  //sp 61
  sp_rates[60] -= fwd_rates[1473];
  //sp 48
  sp_rates[47] -= fwd_rates[1473];

  //rxn 1474
  //sp 49
  sp_rates[48] += fwd_rates[1474];
  //sp 59
  sp_rates[58] += fwd_rates[1474];
  //sp 60
  sp_rates[59] -= fwd_rates[1474];
  //sp 48
  sp_rates[47] -= fwd_rates[1474];

  //rxn 1475
  //sp 49
  sp_rates[48] += fwd_rates[1475];
  //sp 62
  sp_rates[61] += fwd_rates[1475];
  //sp 64
  sp_rates[63] -= fwd_rates[1475];
  //sp 48
  sp_rates[47] -= fwd_rates[1475];

  //rxn 1476
  //sp 65
  sp_rates[64] -= fwd_rates[1476];
  //sp 49
  sp_rates[48] += fwd_rates[1476];
  //sp 62
  sp_rates[61] += fwd_rates[1476];
  //sp 48
  sp_rates[47] -= fwd_rates[1476];

  //rxn 1477
  //sp 49
  sp_rates[48] += fwd_rates[1477];
  //sp 62
  sp_rates[61] += fwd_rates[1477];
  //sp 63
  sp_rates[62] -= fwd_rates[1477];
  //sp 48
  sp_rates[47] -= fwd_rates[1477];

  //rxn 1478
  //sp 49
  sp_rates[48] += fwd_rates[1478];
  //sp 102
  sp_rates[101] += fwd_rates[1478];
  //sp 103
  sp_rates[102] -= fwd_rates[1478];
  //sp 48
  sp_rates[47] -= fwd_rates[1478];

  //rxn 1479
  //sp 49
  sp_rates[48] += fwd_rates[1479];
  //sp 100
  sp_rates[99] -= fwd_rates[1479];
  //sp 97
  sp_rates[96] += fwd_rates[1479];
  //sp 48
  sp_rates[47] -= fwd_rates[1479];

  //rxn 1480
  //sp 49
  sp_rates[48] += fwd_rates[1480];
  //sp 98
  sp_rates[97] += fwd_rates[1480];
  //sp 99
  sp_rates[98] -= fwd_rates[1480];
  //sp 48
  sp_rates[47] -= fwd_rates[1480];

  //rxn 1481
  //sp 49
  sp_rates[48] += fwd_rates[1481];
  //sp 94
  sp_rates[93] += fwd_rates[1481];
  //sp 95
  sp_rates[94] -= fwd_rates[1481];
  //sp 48
  sp_rates[47] -= fwd_rates[1481];

  //rxn 1482
  //sp 105
  sp_rates[104] -= fwd_rates[1482];
  //sp 49
  sp_rates[48] += fwd_rates[1482];
  //sp 104
  sp_rates[103] += fwd_rates[1482];
  //sp 48
  sp_rates[47] -= fwd_rates[1482];

  //rxn 1483
  //sp 49
  sp_rates[48] += fwd_rates[1483];
  //sp 107
  sp_rates[106] += fwd_rates[1483];
  //sp 109
  sp_rates[108] -= fwd_rates[1483];
  //sp 48
  sp_rates[47] -= fwd_rates[1483];

  //rxn 1484
  //sp 49
  sp_rates[48] += fwd_rates[1484];
  //sp 110
  sp_rates[109] += fwd_rates[1484];
  //sp 111
  sp_rates[110] -= fwd_rates[1484];
  //sp 48
  sp_rates[47] -= fwd_rates[1484];

  //rxn 1485
  //sp 49
  sp_rates[48] += fwd_rates[1485];
  //sp 8
  sp_rates[7] += fwd_rates[1485];
  //sp 6
  sp_rates[5] -= fwd_rates[1485];
  //sp 48
  sp_rates[47] -= fwd_rates[1485];

  //rxn 1486
  //sp 49
  sp_rates[48] += fwd_rates[1486];
  //sp 10
  sp_rates[9] += fwd_rates[1486];
  //sp 5
  sp_rates[4] -= fwd_rates[1486];
  //sp 48
  sp_rates[47] -= fwd_rates[1486];

  //rxn 1487
  //sp 49
  sp_rates[48] += fwd_rates[1487];
  //sp 48
  sp_rates[47] -= fwd_rates[1487];
  //sp 7
  sp_rates[6] += fwd_rates[1487];
  //sp 8
  sp_rates[7] -= fwd_rates[1487];

  //rxn 1488
  //sp 9
  sp_rates[8] += fwd_rates[1488];
  //sp 10
  sp_rates[9] -= fwd_rates[1488];
  //sp 49
  sp_rates[48] += fwd_rates[1488];
  //sp 48
  sp_rates[47] -= fwd_rates[1488];

  //rxn 1489
  //sp 49
  sp_rates[48] += fwd_rates[1489];
  //sp 22
  sp_rates[21] += fwd_rates[1489];
  //sp 23
  sp_rates[22] -= fwd_rates[1489];
  //sp 48
  sp_rates[47] -= fwd_rates[1489];

  //rxn 1490
  //sp 49
  sp_rates[48] += fwd_rates[1490];
  //sp 22
  sp_rates[21] += fwd_rates[1490];
  //sp 48
  sp_rates[47] -= fwd_rates[1490];
  //sp 24
  sp_rates[23] -= fwd_rates[1490];

  //rxn 1491
  //sp 25
  sp_rates[24] += fwd_rates[1491];
  //sp 26
  sp_rates[25] -= fwd_rates[1491];
  //sp 49
  sp_rates[48] += fwd_rates[1491];
  //sp 48
  sp_rates[47] -= fwd_rates[1491];

  //rxn 1492
  //sp 49
  sp_rates[48] += fwd_rates[1492];
  //sp 45
  sp_rates[44] += fwd_rates[1492];
  //sp 46
  sp_rates[45] -= fwd_rates[1492];
  //sp 48
  sp_rates[47] -= fwd_rates[1492];

  //rxn 1493
  //sp 49
  sp_rates[48] += fwd_rates[1493];
  //sp 45
  sp_rates[44] += fwd_rates[1493];
  //sp 47
  sp_rates[46] -= fwd_rates[1493];
  //sp 48
  sp_rates[47] -= fwd_rates[1493];

  //rxn 1494
  //sp 49
  sp_rates[48] += fwd_rates[1494];
  //sp 79
  sp_rates[78] += fwd_rates[1494];
  //sp 80
  sp_rates[79] -= fwd_rates[1494];
  //sp 48
  sp_rates[47] -= fwd_rates[1494];

  //rxn 1495
  //sp 3
  sp_rates[2] += fwd_rates[1495];
  //sp 4
  sp_rates[3] -= fwd_rates[1495];
  //sp 11
  sp_rates[10] += fwd_rates[1495];
  //sp 53
  sp_rates[52] -= fwd_rates[1495];
  //sp 24
  sp_rates[23] += fwd_rates[1495];

  //rxn 1496
  //sp 11
  sp_rates[10] += fwd_rates[1496];
  //sp 14
  sp_rates[13] += fwd_rates[1496];
  //sp 15
  sp_rates[14] -= fwd_rates[1496];
  //sp 53
  sp_rates[52] -= fwd_rates[1496];
  //sp 24
  sp_rates[23] += fwd_rates[1496];

  //rxn 1497
  //sp 34
  sp_rates[33] += fwd_rates[1497];
  //sp 35
  sp_rates[34] -= fwd_rates[1497];
  //sp 11
  sp_rates[10] += fwd_rates[1497];
  //sp 53
  sp_rates[52] -= fwd_rates[1497];
  //sp 24
  sp_rates[23] += fwd_rates[1497];

  //rxn 1498
  //sp 11
  sp_rates[10] += fwd_rates[1498];
  //sp 53
  sp_rates[52] -= fwd_rates[1498];
  //sp 24
  sp_rates[23] += fwd_rates[1498];
  //sp 30
  sp_rates[29] += fwd_rates[1498];
  //sp 31
  sp_rates[30] -= fwd_rates[1498];

  //rxn 1499
  //sp 11
  sp_rates[10] += fwd_rates[1499];
  //sp 53
  sp_rates[52] -= fwd_rates[1499];
  //sp 24
  sp_rates[23] += fwd_rates[1499];
  //sp 59
  sp_rates[58] += fwd_rates[1499];
  //sp 61
  sp_rates[60] -= fwd_rates[1499];

  //rxn 1500
  //sp 11
  sp_rates[10] += fwd_rates[1500];
  //sp 53
  sp_rates[52] -= fwd_rates[1500];
  //sp 24
  sp_rates[23] += fwd_rates[1500];
  //sp 59
  sp_rates[58] += fwd_rates[1500];
  //sp 60
  sp_rates[59] -= fwd_rates[1500];

  //rxn 1501
  //sp 11
  sp_rates[10] += fwd_rates[1501];
  //sp 53
  sp_rates[52] -= fwd_rates[1501];
  //sp 24
  sp_rates[23] += fwd_rates[1501];
  //sp 62
  sp_rates[61] += fwd_rates[1501];
  //sp 64
  sp_rates[63] -= fwd_rates[1501];

  //rxn 1502
  //sp 65
  sp_rates[64] -= fwd_rates[1502];
  //sp 11
  sp_rates[10] += fwd_rates[1502];
  //sp 53
  sp_rates[52] -= fwd_rates[1502];
  //sp 24
  sp_rates[23] += fwd_rates[1502];
  //sp 62
  sp_rates[61] += fwd_rates[1502];

  //rxn 1503
  //sp 11
  sp_rates[10] += fwd_rates[1503];
  //sp 53
  sp_rates[52] -= fwd_rates[1503];
  //sp 24
  sp_rates[23] += fwd_rates[1503];
  //sp 62
  sp_rates[61] += fwd_rates[1503];
  //sp 63
  sp_rates[62] -= fwd_rates[1503];

  //rxn 1504
  //sp 102
  sp_rates[101] += fwd_rates[1504];
  //sp 103
  sp_rates[102] -= fwd_rates[1504];
  //sp 11
  sp_rates[10] += fwd_rates[1504];
  //sp 53
  sp_rates[52] -= fwd_rates[1504];
  //sp 24
  sp_rates[23] += fwd_rates[1504];

  //rxn 1505
  //sp 97
  sp_rates[96] += fwd_rates[1505];
  //sp 100
  sp_rates[99] -= fwd_rates[1505];
  //sp 11
  sp_rates[10] += fwd_rates[1505];
  //sp 53
  sp_rates[52] -= fwd_rates[1505];
  //sp 24
  sp_rates[23] += fwd_rates[1505];

  //rxn 1506
  //sp 98
  sp_rates[97] += fwd_rates[1506];
  //sp 99
  sp_rates[98] -= fwd_rates[1506];
  //sp 11
  sp_rates[10] += fwd_rates[1506];
  //sp 53
  sp_rates[52] -= fwd_rates[1506];
  //sp 24
  sp_rates[23] += fwd_rates[1506];

  //rxn 1507
  //sp 11
  sp_rates[10] += fwd_rates[1507];
  //sp 53
  sp_rates[52] -= fwd_rates[1507];
  //sp 24
  sp_rates[23] += fwd_rates[1507];
  //sp 94
  sp_rates[93] += fwd_rates[1507];
  //sp 95
  sp_rates[94] -= fwd_rates[1507];

  //rxn 1508
  //sp 104
  sp_rates[103] += fwd_rates[1508];
  //sp 105
  sp_rates[104] -= fwd_rates[1508];
  //sp 11
  sp_rates[10] += fwd_rates[1508];
  //sp 53
  sp_rates[52] -= fwd_rates[1508];
  //sp 24
  sp_rates[23] += fwd_rates[1508];

  //rxn 1509
  //sp 11
  sp_rates[10] += fwd_rates[1509];
  //sp 107
  sp_rates[106] += fwd_rates[1509];
  //sp 109
  sp_rates[108] -= fwd_rates[1509];
  //sp 53
  sp_rates[52] -= fwd_rates[1509];
  //sp 24
  sp_rates[23] += fwd_rates[1509];

  //rxn 1510
  //sp 11
  sp_rates[10] += fwd_rates[1510];
  //sp 110
  sp_rates[109] += fwd_rates[1510];
  //sp 111
  sp_rates[110] -= fwd_rates[1510];
  //sp 53
  sp_rates[52] -= fwd_rates[1510];
  //sp 24
  sp_rates[23] += fwd_rates[1510];

  //rxn 1511
  //sp 6
  sp_rates[5] -= fwd_rates[1511];
  //sp 8
  sp_rates[7] += fwd_rates[1511];
  //sp 11
  sp_rates[10] += fwd_rates[1511];
  //sp 53
  sp_rates[52] -= fwd_rates[1511];
  //sp 24
  sp_rates[23] += fwd_rates[1511];

  //rxn 1512
  //sp 5
  sp_rates[4] -= fwd_rates[1512];
  //sp 10
  sp_rates[9] += fwd_rates[1512];
  //sp 11
  sp_rates[10] += fwd_rates[1512];
  //sp 53
  sp_rates[52] -= fwd_rates[1512];
  //sp 24
  sp_rates[23] += fwd_rates[1512];

  //rxn 1513
  //sp 7
  sp_rates[6] += fwd_rates[1513];
  //sp 8
  sp_rates[7] -= fwd_rates[1513];
  //sp 11
  sp_rates[10] += fwd_rates[1513];
  //sp 53
  sp_rates[52] -= fwd_rates[1513];
  //sp 24
  sp_rates[23] += fwd_rates[1513];

  //rxn 1514
  //sp 9
  sp_rates[8] += fwd_rates[1514];
  //sp 10
  sp_rates[9] -= fwd_rates[1514];
  //sp 11
  sp_rates[10] += fwd_rates[1514];
  //sp 53
  sp_rates[52] -= fwd_rates[1514];
  //sp 24
  sp_rates[23] += fwd_rates[1514];

  //rxn 1515
  //sp 11
  sp_rates[10] += fwd_rates[1515];
  //sp 53
  sp_rates[52] -= fwd_rates[1515];
  //sp 22
  sp_rates[21] += fwd_rates[1515];
  //sp 23
  sp_rates[22] -= fwd_rates[1515];
  //sp 24
  sp_rates[23] += fwd_rates[1515];

  //rxn 1516
  //sp 11
  sp_rates[10] += fwd_rates[1516];
  //sp 53
  sp_rates[52] -= fwd_rates[1516];
  //sp 22
  sp_rates[21] += fwd_rates[1516];

  //rxn 1517
  //sp 11
  sp_rates[10] += fwd_rates[1517];
  //sp 53
  sp_rates[52] -= fwd_rates[1517];
  //sp 24
  sp_rates[23] += fwd_rates[1517];
  //sp 25
  sp_rates[24] += fwd_rates[1517];
  //sp 26
  sp_rates[25] -= fwd_rates[1517];

  //rxn 1518
  //sp 11
  sp_rates[10] += fwd_rates[1518];
  //sp 45
  sp_rates[44] += fwd_rates[1518];
  //sp 46
  sp_rates[45] -= fwd_rates[1518];
  //sp 53
  sp_rates[52] -= fwd_rates[1518];
  //sp 24
  sp_rates[23] += fwd_rates[1518];

  //rxn 1519
  //sp 11
  sp_rates[10] += fwd_rates[1519];
  //sp 45
  sp_rates[44] += fwd_rates[1519];
  //sp 47
  sp_rates[46] -= fwd_rates[1519];
  //sp 53
  sp_rates[52] -= fwd_rates[1519];
  //sp 24
  sp_rates[23] += fwd_rates[1519];

  //rxn 1520
  //sp 11
  sp_rates[10] += fwd_rates[1520];
  //sp 79
  sp_rates[78] += fwd_rates[1520];
  //sp 80
  sp_rates[79] -= fwd_rates[1520];
  //sp 53
  sp_rates[52] -= fwd_rates[1520];
  //sp 24
  sp_rates[23] += fwd_rates[1520];

  //rxn 1521
  //sp 11
  sp_rates[10] += fwd_rates[1521];
  //sp 48
  sp_rates[47] += fwd_rates[1521];
  //sp 49
  sp_rates[48] -= fwd_rates[1521];
  //sp 53
  sp_rates[52] -= fwd_rates[1521];
  //sp 24
  sp_rates[23] += fwd_rates[1521];

  //rxn 1522
  //sp 3
  sp_rates[2] += fwd_rates[1522];
  //sp 53
  sp_rates[52] -= fwd_rates[1522];
  //sp 54
  sp_rates[53] += fwd_rates[1522];

  //rxn 1523
  //sp 4
  sp_rates[3] += fwd_rates[1523];
  //sp 14
  sp_rates[13] += fwd_rates[1523];
  //sp 15
  sp_rates[14] -= fwd_rates[1523];
  //sp 53
  sp_rates[52] -= fwd_rates[1523];
  //sp 54
  sp_rates[53] += fwd_rates[1523];

  //rxn 1524
  //sp 34
  sp_rates[33] += fwd_rates[1524];
  //sp 35
  sp_rates[34] -= fwd_rates[1524];
  //sp 4
  sp_rates[3] += fwd_rates[1524];
  //sp 53
  sp_rates[52] -= fwd_rates[1524];
  //sp 54
  sp_rates[53] += fwd_rates[1524];

  //rxn 1525
  //sp 4
  sp_rates[3] += fwd_rates[1525];
  //sp 53
  sp_rates[52] -= fwd_rates[1525];
  //sp 54
  sp_rates[53] += fwd_rates[1525];
  //sp 30
  sp_rates[29] += fwd_rates[1525];
  //sp 31
  sp_rates[30] -= fwd_rates[1525];

  //rxn 1526
  //sp 4
  sp_rates[3] += fwd_rates[1526];
  //sp 53
  sp_rates[52] -= fwd_rates[1526];
  //sp 54
  sp_rates[53] += fwd_rates[1526];
  //sp 59
  sp_rates[58] += fwd_rates[1526];
  //sp 61
  sp_rates[60] -= fwd_rates[1526];

  //rxn 1527
  //sp 4
  sp_rates[3] += fwd_rates[1527];
  //sp 53
  sp_rates[52] -= fwd_rates[1527];
  //sp 54
  sp_rates[53] += fwd_rates[1527];
  //sp 59
  sp_rates[58] += fwd_rates[1527];
  //sp 60
  sp_rates[59] -= fwd_rates[1527];

  //rxn 1528
  //sp 4
  sp_rates[3] += fwd_rates[1528];
  //sp 53
  sp_rates[52] -= fwd_rates[1528];
  //sp 54
  sp_rates[53] += fwd_rates[1528];
  //sp 62
  sp_rates[61] += fwd_rates[1528];
  //sp 64
  sp_rates[63] -= fwd_rates[1528];

  //rxn 1529
  //sp 65
  sp_rates[64] -= fwd_rates[1529];
  //sp 4
  sp_rates[3] += fwd_rates[1529];
  //sp 53
  sp_rates[52] -= fwd_rates[1529];
  //sp 54
  sp_rates[53] += fwd_rates[1529];
  //sp 62
  sp_rates[61] += fwd_rates[1529];

  //rxn 1530
  //sp 4
  sp_rates[3] += fwd_rates[1530];
  //sp 53
  sp_rates[52] -= fwd_rates[1530];
  //sp 54
  sp_rates[53] += fwd_rates[1530];
  //sp 62
  sp_rates[61] += fwd_rates[1530];
  //sp 63
  sp_rates[62] -= fwd_rates[1530];

  //rxn 1531
  //sp 4
  sp_rates[3] += fwd_rates[1531];
  //sp 102
  sp_rates[101] += fwd_rates[1531];
  //sp 103
  sp_rates[102] -= fwd_rates[1531];
  //sp 53
  sp_rates[52] -= fwd_rates[1531];
  //sp 54
  sp_rates[53] += fwd_rates[1531];

  //rxn 1532
  //sp 97
  sp_rates[96] += fwd_rates[1532];
  //sp 4
  sp_rates[3] += fwd_rates[1532];
  //sp 100
  sp_rates[99] -= fwd_rates[1532];
  //sp 53
  sp_rates[52] -= fwd_rates[1532];
  //sp 54
  sp_rates[53] += fwd_rates[1532];

  //rxn 1533
  //sp 98
  sp_rates[97] += fwd_rates[1533];
  //sp 99
  sp_rates[98] -= fwd_rates[1533];
  //sp 4
  sp_rates[3] += fwd_rates[1533];
  //sp 53
  sp_rates[52] -= fwd_rates[1533];
  //sp 54
  sp_rates[53] += fwd_rates[1533];

  //rxn 1534
  //sp 4
  sp_rates[3] += fwd_rates[1534];
  //sp 53
  sp_rates[52] -= fwd_rates[1534];
  //sp 54
  sp_rates[53] += fwd_rates[1534];
  //sp 94
  sp_rates[93] += fwd_rates[1534];
  //sp 95
  sp_rates[94] -= fwd_rates[1534];

  //rxn 1535
  //sp 4
  sp_rates[3] += fwd_rates[1535];
  //sp 104
  sp_rates[103] += fwd_rates[1535];
  //sp 105
  sp_rates[104] -= fwd_rates[1535];
  //sp 53
  sp_rates[52] -= fwd_rates[1535];
  //sp 54
  sp_rates[53] += fwd_rates[1535];

  //rxn 1536
  //sp 4
  sp_rates[3] += fwd_rates[1536];
  //sp 107
  sp_rates[106] += fwd_rates[1536];
  //sp 109
  sp_rates[108] -= fwd_rates[1536];
  //sp 53
  sp_rates[52] -= fwd_rates[1536];
  //sp 54
  sp_rates[53] += fwd_rates[1536];

  //rxn 1537
  //sp 4
  sp_rates[3] += fwd_rates[1537];
  //sp 110
  sp_rates[109] += fwd_rates[1537];
  //sp 111
  sp_rates[110] -= fwd_rates[1537];
  //sp 53
  sp_rates[52] -= fwd_rates[1537];
  //sp 54
  sp_rates[53] += fwd_rates[1537];

  //rxn 1538
  //sp 4
  sp_rates[3] += fwd_rates[1538];
  //sp 6
  sp_rates[5] -= fwd_rates[1538];
  //sp 8
  sp_rates[7] += fwd_rates[1538];
  //sp 53
  sp_rates[52] -= fwd_rates[1538];
  //sp 54
  sp_rates[53] += fwd_rates[1538];

  //rxn 1539
  //sp 4
  sp_rates[3] += fwd_rates[1539];
  //sp 5
  sp_rates[4] -= fwd_rates[1539];
  //sp 10
  sp_rates[9] += fwd_rates[1539];
  //sp 53
  sp_rates[52] -= fwd_rates[1539];
  //sp 54
  sp_rates[53] += fwd_rates[1539];

  //rxn 1540
  //sp 4
  sp_rates[3] += fwd_rates[1540];
  //sp 7
  sp_rates[6] += fwd_rates[1540];
  //sp 8
  sp_rates[7] -= fwd_rates[1540];
  //sp 53
  sp_rates[52] -= fwd_rates[1540];
  //sp 54
  sp_rates[53] += fwd_rates[1540];

  //rxn 1541
  //sp 4
  sp_rates[3] += fwd_rates[1541];
  //sp 9
  sp_rates[8] += fwd_rates[1541];
  //sp 10
  sp_rates[9] -= fwd_rates[1541];
  //sp 53
  sp_rates[52] -= fwd_rates[1541];
  //sp 54
  sp_rates[53] += fwd_rates[1541];

  //rxn 1542
  //sp 4
  sp_rates[3] += fwd_rates[1542];
  //sp 53
  sp_rates[52] -= fwd_rates[1542];
  //sp 54
  sp_rates[53] += fwd_rates[1542];
  //sp 22
  sp_rates[21] += fwd_rates[1542];
  //sp 23
  sp_rates[22] -= fwd_rates[1542];

  //rxn 1543
  //sp 4
  sp_rates[3] += fwd_rates[1543];
  //sp 53
  sp_rates[52] -= fwd_rates[1543];
  //sp 54
  sp_rates[53] += fwd_rates[1543];
  //sp 22
  sp_rates[21] += fwd_rates[1543];
  //sp 24
  sp_rates[23] -= fwd_rates[1543];

  //rxn 1544
  //sp 4
  sp_rates[3] += fwd_rates[1544];
  //sp 53
  sp_rates[52] -= fwd_rates[1544];
  //sp 54
  sp_rates[53] += fwd_rates[1544];
  //sp 25
  sp_rates[24] += fwd_rates[1544];
  //sp 26
  sp_rates[25] -= fwd_rates[1544];

  //rxn 1545
  //sp 4
  sp_rates[3] += fwd_rates[1545];
  //sp 45
  sp_rates[44] += fwd_rates[1545];
  //sp 46
  sp_rates[45] -= fwd_rates[1545];
  //sp 53
  sp_rates[52] -= fwd_rates[1545];
  //sp 54
  sp_rates[53] += fwd_rates[1545];

  //rxn 1546
  //sp 4
  sp_rates[3] += fwd_rates[1546];
  //sp 45
  sp_rates[44] += fwd_rates[1546];
  //sp 47
  sp_rates[46] -= fwd_rates[1546];
  //sp 53
  sp_rates[52] -= fwd_rates[1546];
  //sp 54
  sp_rates[53] += fwd_rates[1546];

  //rxn 1547
  //sp 4
  sp_rates[3] += fwd_rates[1547];
  //sp 79
  sp_rates[78] += fwd_rates[1547];
  //sp 80
  sp_rates[79] -= fwd_rates[1547];
  //sp 53
  sp_rates[52] -= fwd_rates[1547];
  //sp 54
  sp_rates[53] += fwd_rates[1547];

  //rxn 1548
  //sp 4
  sp_rates[3] += fwd_rates[1548];
  //sp 48
  sp_rates[47] += fwd_rates[1548];
  //sp 49
  sp_rates[48] -= fwd_rates[1548];
  //sp 53
  sp_rates[52] -= fwd_rates[1548];
  //sp 54
  sp_rates[53] += fwd_rates[1548];

  //rxn 1549
  //sp 3
  sp_rates[2] += fwd_rates[1549];
  //sp 4
  sp_rates[3] -= fwd_rates[1549];
  //sp 11
  sp_rates[10] += fwd_rates[1549];
  //sp 54
  sp_rates[53] -= fwd_rates[1549];
  //sp 26
  sp_rates[25] += fwd_rates[1549];

  //rxn 1550
  //sp 11
  sp_rates[10] += fwd_rates[1550];
  //sp 14
  sp_rates[13] += fwd_rates[1550];
  //sp 15
  sp_rates[14] -= fwd_rates[1550];
  //sp 54
  sp_rates[53] -= fwd_rates[1550];
  //sp 26
  sp_rates[25] += fwd_rates[1550];

  //rxn 1551
  //sp 34
  sp_rates[33] += fwd_rates[1551];
  //sp 35
  sp_rates[34] -= fwd_rates[1551];
  //sp 11
  sp_rates[10] += fwd_rates[1551];
  //sp 54
  sp_rates[53] -= fwd_rates[1551];
  //sp 26
  sp_rates[25] += fwd_rates[1551];

  //rxn 1552
  //sp 11
  sp_rates[10] += fwd_rates[1552];
  //sp 54
  sp_rates[53] -= fwd_rates[1552];
  //sp 26
  sp_rates[25] += fwd_rates[1552];
  //sp 30
  sp_rates[29] += fwd_rates[1552];
  //sp 31
  sp_rates[30] -= fwd_rates[1552];

  //rxn 1553
  //sp 11
  sp_rates[10] += fwd_rates[1553];
  //sp 54
  sp_rates[53] -= fwd_rates[1553];
  //sp 26
  sp_rates[25] += fwd_rates[1553];
  //sp 59
  sp_rates[58] += fwd_rates[1553];
  //sp 61
  sp_rates[60] -= fwd_rates[1553];

  //rxn 1554
  //sp 11
  sp_rates[10] += fwd_rates[1554];
  //sp 54
  sp_rates[53] -= fwd_rates[1554];
  //sp 26
  sp_rates[25] += fwd_rates[1554];
  //sp 59
  sp_rates[58] += fwd_rates[1554];
  //sp 60
  sp_rates[59] -= fwd_rates[1554];

  //rxn 1555
  //sp 11
  sp_rates[10] += fwd_rates[1555];
  //sp 54
  sp_rates[53] -= fwd_rates[1555];
  //sp 26
  sp_rates[25] += fwd_rates[1555];
  //sp 62
  sp_rates[61] += fwd_rates[1555];
  //sp 64
  sp_rates[63] -= fwd_rates[1555];

  //rxn 1556
  //sp 65
  sp_rates[64] -= fwd_rates[1556];
  //sp 11
  sp_rates[10] += fwd_rates[1556];
  //sp 54
  sp_rates[53] -= fwd_rates[1556];
  //sp 26
  sp_rates[25] += fwd_rates[1556];
  //sp 62
  sp_rates[61] += fwd_rates[1556];

  //rxn 1557
  //sp 11
  sp_rates[10] += fwd_rates[1557];
  //sp 54
  sp_rates[53] -= fwd_rates[1557];
  //sp 26
  sp_rates[25] += fwd_rates[1557];
  //sp 62
  sp_rates[61] += fwd_rates[1557];
  //sp 63
  sp_rates[62] -= fwd_rates[1557];

  //rxn 1558
  //sp 102
  sp_rates[101] += fwd_rates[1558];
  //sp 103
  sp_rates[102] -= fwd_rates[1558];
  //sp 11
  sp_rates[10] += fwd_rates[1558];
  //sp 54
  sp_rates[53] -= fwd_rates[1558];
  //sp 26
  sp_rates[25] += fwd_rates[1558];

  //rxn 1559
  //sp 97
  sp_rates[96] += fwd_rates[1559];
  //sp 100
  sp_rates[99] -= fwd_rates[1559];
  //sp 11
  sp_rates[10] += fwd_rates[1559];
  //sp 54
  sp_rates[53] -= fwd_rates[1559];
  //sp 26
  sp_rates[25] += fwd_rates[1559];

  //rxn 1560
  //sp 98
  sp_rates[97] += fwd_rates[1560];
  //sp 99
  sp_rates[98] -= fwd_rates[1560];
  //sp 11
  sp_rates[10] += fwd_rates[1560];
  //sp 54
  sp_rates[53] -= fwd_rates[1560];
  //sp 26
  sp_rates[25] += fwd_rates[1560];

  //rxn 1561
  //sp 11
  sp_rates[10] += fwd_rates[1561];
  //sp 54
  sp_rates[53] -= fwd_rates[1561];
  //sp 26
  sp_rates[25] += fwd_rates[1561];
  //sp 94
  sp_rates[93] += fwd_rates[1561];
  //sp 95
  sp_rates[94] -= fwd_rates[1561];

  //rxn 1562
  //sp 104
  sp_rates[103] += fwd_rates[1562];
  //sp 105
  sp_rates[104] -= fwd_rates[1562];
  //sp 11
  sp_rates[10] += fwd_rates[1562];
  //sp 54
  sp_rates[53] -= fwd_rates[1562];
  //sp 26
  sp_rates[25] += fwd_rates[1562];

  //rxn 1563
  //sp 11
  sp_rates[10] += fwd_rates[1563];
  //sp 107
  sp_rates[106] += fwd_rates[1563];
  //sp 109
  sp_rates[108] -= fwd_rates[1563];
  //sp 54
  sp_rates[53] -= fwd_rates[1563];
  //sp 26
  sp_rates[25] += fwd_rates[1563];

  //rxn 1564
  //sp 11
  sp_rates[10] += fwd_rates[1564];
  //sp 110
  sp_rates[109] += fwd_rates[1564];
  //sp 111
  sp_rates[110] -= fwd_rates[1564];
  //sp 54
  sp_rates[53] -= fwd_rates[1564];
  //sp 26
  sp_rates[25] += fwd_rates[1564];

  //rxn 1565
  //sp 6
  sp_rates[5] -= fwd_rates[1565];
  //sp 8
  sp_rates[7] += fwd_rates[1565];
  //sp 11
  sp_rates[10] += fwd_rates[1565];
  //sp 54
  sp_rates[53] -= fwd_rates[1565];
  //sp 26
  sp_rates[25] += fwd_rates[1565];

  //rxn 1566
  //sp 5
  sp_rates[4] -= fwd_rates[1566];
  //sp 10
  sp_rates[9] += fwd_rates[1566];
  //sp 11
  sp_rates[10] += fwd_rates[1566];
  //sp 54
  sp_rates[53] -= fwd_rates[1566];
  //sp 26
  sp_rates[25] += fwd_rates[1566];

  //rxn 1567
  //sp 7
  sp_rates[6] += fwd_rates[1567];
  //sp 8
  sp_rates[7] -= fwd_rates[1567];
  //sp 11
  sp_rates[10] += fwd_rates[1567];
  //sp 54
  sp_rates[53] -= fwd_rates[1567];
  //sp 26
  sp_rates[25] += fwd_rates[1567];

  //rxn 1568
  //sp 9
  sp_rates[8] += fwd_rates[1568];
  //sp 10
  sp_rates[9] -= fwd_rates[1568];
  //sp 11
  sp_rates[10] += fwd_rates[1568];
  //sp 54
  sp_rates[53] -= fwd_rates[1568];
  //sp 26
  sp_rates[25] += fwd_rates[1568];

  //rxn 1569
  //sp 11
  sp_rates[10] += fwd_rates[1569];
  //sp 22
  sp_rates[21] += fwd_rates[1569];
  //sp 54
  sp_rates[53] -= fwd_rates[1569];
  //sp 23
  sp_rates[22] -= fwd_rates[1569];
  //sp 26
  sp_rates[25] += fwd_rates[1569];

  //rxn 1570
  //sp 11
  sp_rates[10] += fwd_rates[1570];
  //sp 22
  sp_rates[21] += fwd_rates[1570];
  //sp 54
  sp_rates[53] -= fwd_rates[1570];
  //sp 24
  sp_rates[23] -= fwd_rates[1570];
  //sp 26
  sp_rates[25] += fwd_rates[1570];

  //rxn 1571
  //sp 25
  sp_rates[24] += fwd_rates[1571];
  //sp 11
  sp_rates[10] += fwd_rates[1571];
  //sp 54
  sp_rates[53] -= fwd_rates[1571];

  //rxn 1572
  //sp 11
  sp_rates[10] += fwd_rates[1572];
  //sp 45
  sp_rates[44] += fwd_rates[1572];
  //sp 46
  sp_rates[45] -= fwd_rates[1572];
  //sp 54
  sp_rates[53] -= fwd_rates[1572];
  //sp 26
  sp_rates[25] += fwd_rates[1572];

  //rxn 1573
  //sp 11
  sp_rates[10] += fwd_rates[1573];
  //sp 45
  sp_rates[44] += fwd_rates[1573];
  //sp 47
  sp_rates[46] -= fwd_rates[1573];
  //sp 54
  sp_rates[53] -= fwd_rates[1573];
  //sp 26
  sp_rates[25] += fwd_rates[1573];

  //rxn 1574
  //sp 11
  sp_rates[10] += fwd_rates[1574];
  //sp 79
  sp_rates[78] += fwd_rates[1574];
  //sp 80
  sp_rates[79] -= fwd_rates[1574];
  //sp 54
  sp_rates[53] -= fwd_rates[1574];
  //sp 26
  sp_rates[25] += fwd_rates[1574];

  //rxn 1575
  //sp 11
  sp_rates[10] += fwd_rates[1575];
  //sp 48
  sp_rates[47] += fwd_rates[1575];
  //sp 49
  sp_rates[48] -= fwd_rates[1575];
  //sp 54
  sp_rates[53] -= fwd_rates[1575];
  //sp 26
  sp_rates[25] += fwd_rates[1575];

  //rxn 1576
  //sp 3
  sp_rates[2] += fwd_rates[1576];
  //sp 4
  sp_rates[3] -= fwd_rates[1576];
  //sp 8
  sp_rates[7] += fwd_rates[1576];
  //sp 11
  sp_rates[10] += fwd_rates[1576];
  //sp 50
  sp_rates[49] -= fwd_rates[1576];

  //rxn 1577
  //sp 8
  sp_rates[7] += fwd_rates[1577];
  //sp 11
  sp_rates[10] += fwd_rates[1577];
  //sp 14
  sp_rates[13] += fwd_rates[1577];
  //sp 15
  sp_rates[14] -= fwd_rates[1577];
  //sp 50
  sp_rates[49] -= fwd_rates[1577];

  //rxn 1578
  //sp 34
  sp_rates[33] += fwd_rates[1578];
  //sp 35
  sp_rates[34] -= fwd_rates[1578];
  //sp 8
  sp_rates[7] += fwd_rates[1578];
  //sp 11
  sp_rates[10] += fwd_rates[1578];
  //sp 50
  sp_rates[49] -= fwd_rates[1578];

  //rxn 1579
  //sp 8
  sp_rates[7] += fwd_rates[1579];
  //sp 11
  sp_rates[10] += fwd_rates[1579];
  //sp 50
  sp_rates[49] -= fwd_rates[1579];
  //sp 30
  sp_rates[29] += fwd_rates[1579];
  //sp 31
  sp_rates[30] -= fwd_rates[1579];

  //rxn 1580
  //sp 8
  sp_rates[7] += fwd_rates[1580];
  //sp 11
  sp_rates[10] += fwd_rates[1580];
  //sp 50
  sp_rates[49] -= fwd_rates[1580];
  //sp 59
  sp_rates[58] += fwd_rates[1580];
  //sp 61
  sp_rates[60] -= fwd_rates[1580];

  //rxn 1581
  //sp 8
  sp_rates[7] += fwd_rates[1581];
  //sp 11
  sp_rates[10] += fwd_rates[1581];
  //sp 50
  sp_rates[49] -= fwd_rates[1581];
  //sp 59
  sp_rates[58] += fwd_rates[1581];
  //sp 60
  sp_rates[59] -= fwd_rates[1581];

  //rxn 1582
  //sp 8
  sp_rates[7] += fwd_rates[1582];
  //sp 11
  sp_rates[10] += fwd_rates[1582];
  //sp 50
  sp_rates[49] -= fwd_rates[1582];
  //sp 62
  sp_rates[61] += fwd_rates[1582];
  //sp 64
  sp_rates[63] -= fwd_rates[1582];

  //rxn 1583
  //sp 65
  sp_rates[64] -= fwd_rates[1583];
  //sp 8
  sp_rates[7] += fwd_rates[1583];
  //sp 11
  sp_rates[10] += fwd_rates[1583];
  //sp 50
  sp_rates[49] -= fwd_rates[1583];
  //sp 62
  sp_rates[61] += fwd_rates[1583];

  //rxn 1584
  //sp 8
  sp_rates[7] += fwd_rates[1584];
  //sp 11
  sp_rates[10] += fwd_rates[1584];
  //sp 50
  sp_rates[49] -= fwd_rates[1584];
  //sp 62
  sp_rates[61] += fwd_rates[1584];
  //sp 63
  sp_rates[62] -= fwd_rates[1584];

  //rxn 1585
  //sp 102
  sp_rates[101] += fwd_rates[1585];
  //sp 103
  sp_rates[102] -= fwd_rates[1585];
  //sp 8
  sp_rates[7] += fwd_rates[1585];
  //sp 11
  sp_rates[10] += fwd_rates[1585];
  //sp 50
  sp_rates[49] -= fwd_rates[1585];

  //rxn 1586
  //sp 97
  sp_rates[96] += fwd_rates[1586];
  //sp 100
  sp_rates[99] -= fwd_rates[1586];
  //sp 8
  sp_rates[7] += fwd_rates[1586];
  //sp 11
  sp_rates[10] += fwd_rates[1586];
  //sp 50
  sp_rates[49] -= fwd_rates[1586];

  //rxn 1587
  //sp 98
  sp_rates[97] += fwd_rates[1587];
  //sp 99
  sp_rates[98] -= fwd_rates[1587];
  //sp 8
  sp_rates[7] += fwd_rates[1587];
  //sp 11
  sp_rates[10] += fwd_rates[1587];
  //sp 50
  sp_rates[49] -= fwd_rates[1587];

  //rxn 1588
  //sp 8
  sp_rates[7] += fwd_rates[1588];
  //sp 11
  sp_rates[10] += fwd_rates[1588];
  //sp 50
  sp_rates[49] -= fwd_rates[1588];
  //sp 94
  sp_rates[93] += fwd_rates[1588];
  //sp 95
  sp_rates[94] -= fwd_rates[1588];

  //rxn 1589
  //sp 104
  sp_rates[103] += fwd_rates[1589];
  //sp 105
  sp_rates[104] -= fwd_rates[1589];
  //sp 8
  sp_rates[7] += fwd_rates[1589];
  //sp 11
  sp_rates[10] += fwd_rates[1589];
  //sp 50
  sp_rates[49] -= fwd_rates[1589];

  //rxn 1590
  //sp 8
  sp_rates[7] += fwd_rates[1590];
  //sp 11
  sp_rates[10] += fwd_rates[1590];
  //sp 107
  sp_rates[106] += fwd_rates[1590];
  //sp 109
  sp_rates[108] -= fwd_rates[1590];
  //sp 50
  sp_rates[49] -= fwd_rates[1590];

  //rxn 1591
  //sp 8
  sp_rates[7] += fwd_rates[1591];
  //sp 11
  sp_rates[10] += fwd_rates[1591];
  //sp 110
  sp_rates[109] += fwd_rates[1591];
  //sp 111
  sp_rates[110] -= fwd_rates[1591];
  //sp 50
  sp_rates[49] -= fwd_rates[1591];

  //rxn 1592
  //sp 50
  sp_rates[49] -= fwd_rates[1592];
  //sp 11
  sp_rates[10] += fwd_rates[1592];
  //sp 6
  sp_rates[5] -= fwd_rates[1592];
  //sp 8
  sp_rates[7] += 2.0 * fwd_rates[1592];

  //rxn 1593
  //sp 5
  sp_rates[4] -= fwd_rates[1593];
  //sp 8
  sp_rates[7] += fwd_rates[1593];
  //sp 10
  sp_rates[9] += fwd_rates[1593];
  //sp 11
  sp_rates[10] += fwd_rates[1593];
  //sp 50
  sp_rates[49] -= fwd_rates[1593];

  //rxn 1594
  //sp 50
  sp_rates[49] -= fwd_rates[1594];
  //sp 11
  sp_rates[10] += fwd_rates[1594];
  //sp 7
  sp_rates[6] += fwd_rates[1594];

  //rxn 1595
  //sp 8
  sp_rates[7] += fwd_rates[1595];
  //sp 9
  sp_rates[8] += fwd_rates[1595];
  //sp 10
  sp_rates[9] -= fwd_rates[1595];
  //sp 11
  sp_rates[10] += fwd_rates[1595];
  //sp 50
  sp_rates[49] -= fwd_rates[1595];

  //rxn 1596
  //sp 8
  sp_rates[7] += fwd_rates[1596];
  //sp 11
  sp_rates[10] += fwd_rates[1596];
  //sp 50
  sp_rates[49] -= fwd_rates[1596];
  //sp 22
  sp_rates[21] += fwd_rates[1596];
  //sp 23
  sp_rates[22] -= fwd_rates[1596];

  //rxn 1597
  //sp 8
  sp_rates[7] += fwd_rates[1597];
  //sp 11
  sp_rates[10] += fwd_rates[1597];
  //sp 50
  sp_rates[49] -= fwd_rates[1597];
  //sp 22
  sp_rates[21] += fwd_rates[1597];
  //sp 24
  sp_rates[23] -= fwd_rates[1597];

  //rxn 1598
  //sp 8
  sp_rates[7] += fwd_rates[1598];
  //sp 11
  sp_rates[10] += fwd_rates[1598];
  //sp 50
  sp_rates[49] -= fwd_rates[1598];
  //sp 25
  sp_rates[24] += fwd_rates[1598];
  //sp 26
  sp_rates[25] -= fwd_rates[1598];

  //rxn 1599
  //sp 8
  sp_rates[7] += fwd_rates[1599];
  //sp 11
  sp_rates[10] += fwd_rates[1599];
  //sp 45
  sp_rates[44] += fwd_rates[1599];
  //sp 46
  sp_rates[45] -= fwd_rates[1599];
  //sp 50
  sp_rates[49] -= fwd_rates[1599];

  //rxn 1600
  //sp 8
  sp_rates[7] += fwd_rates[1600];
  //sp 11
  sp_rates[10] += fwd_rates[1600];
  //sp 45
  sp_rates[44] += fwd_rates[1600];
  //sp 47
  sp_rates[46] -= fwd_rates[1600];
  //sp 50
  sp_rates[49] -= fwd_rates[1600];

  //rxn 1601
  //sp 8
  sp_rates[7] += fwd_rates[1601];
  //sp 11
  sp_rates[10] += fwd_rates[1601];
  //sp 79
  sp_rates[78] += fwd_rates[1601];
  //sp 80
  sp_rates[79] -= fwd_rates[1601];
  //sp 50
  sp_rates[49] -= fwd_rates[1601];

  //rxn 1602
  //sp 8
  sp_rates[7] += fwd_rates[1602];
  //sp 11
  sp_rates[10] += fwd_rates[1602];
  //sp 48
  sp_rates[47] += fwd_rates[1602];
  //sp 49
  sp_rates[48] -= fwd_rates[1602];
  //sp 50
  sp_rates[49] -= fwd_rates[1602];

  //rxn 1603
  //sp 12
  sp_rates[11] += fwd_rates[1603];
  //sp 50
  sp_rates[49] -= fwd_rates[1603];
  //sp 3
  sp_rates[2] += fwd_rates[1603];

  //rxn 1604
  //sp 4
  sp_rates[3] += fwd_rates[1604];
  //sp 12
  sp_rates[11] += fwd_rates[1604];
  //sp 14
  sp_rates[13] += fwd_rates[1604];
  //sp 15
  sp_rates[14] -= fwd_rates[1604];
  //sp 50
  sp_rates[49] -= fwd_rates[1604];

  //rxn 1605
  //sp 34
  sp_rates[33] += fwd_rates[1605];
  //sp 35
  sp_rates[34] -= fwd_rates[1605];
  //sp 4
  sp_rates[3] += fwd_rates[1605];
  //sp 12
  sp_rates[11] += fwd_rates[1605];
  //sp 50
  sp_rates[49] -= fwd_rates[1605];

  //rxn 1606
  //sp 4
  sp_rates[3] += fwd_rates[1606];
  //sp 12
  sp_rates[11] += fwd_rates[1606];
  //sp 50
  sp_rates[49] -= fwd_rates[1606];
  //sp 30
  sp_rates[29] += fwd_rates[1606];
  //sp 31
  sp_rates[30] -= fwd_rates[1606];

  //rxn 1607
  //sp 4
  sp_rates[3] += fwd_rates[1607];
  //sp 12
  sp_rates[11] += fwd_rates[1607];
  //sp 50
  sp_rates[49] -= fwd_rates[1607];
  //sp 59
  sp_rates[58] += fwd_rates[1607];
  //sp 61
  sp_rates[60] -= fwd_rates[1607];

  //rxn 1608
  //sp 4
  sp_rates[3] += fwd_rates[1608];
  //sp 12
  sp_rates[11] += fwd_rates[1608];
  //sp 50
  sp_rates[49] -= fwd_rates[1608];
  //sp 59
  sp_rates[58] += fwd_rates[1608];
  //sp 60
  sp_rates[59] -= fwd_rates[1608];

  //rxn 1609
  //sp 4
  sp_rates[3] += fwd_rates[1609];
  //sp 12
  sp_rates[11] += fwd_rates[1609];
  //sp 50
  sp_rates[49] -= fwd_rates[1609];
  //sp 62
  sp_rates[61] += fwd_rates[1609];
  //sp 64
  sp_rates[63] -= fwd_rates[1609];

  //rxn 1610
  //sp 65
  sp_rates[64] -= fwd_rates[1610];
  //sp 4
  sp_rates[3] += fwd_rates[1610];
  //sp 12
  sp_rates[11] += fwd_rates[1610];
  //sp 50
  sp_rates[49] -= fwd_rates[1610];
  //sp 62
  sp_rates[61] += fwd_rates[1610];

  //rxn 1611
  //sp 4
  sp_rates[3] += fwd_rates[1611];
  //sp 12
  sp_rates[11] += fwd_rates[1611];
  //sp 50
  sp_rates[49] -= fwd_rates[1611];
  //sp 62
  sp_rates[61] += fwd_rates[1611];
  //sp 63
  sp_rates[62] -= fwd_rates[1611];

  //rxn 1612
  //sp 4
  sp_rates[3] += fwd_rates[1612];
  //sp 102
  sp_rates[101] += fwd_rates[1612];
  //sp 103
  sp_rates[102] -= fwd_rates[1612];
  //sp 12
  sp_rates[11] += fwd_rates[1612];
  //sp 50
  sp_rates[49] -= fwd_rates[1612];

  //rxn 1613
  //sp 97
  sp_rates[96] += fwd_rates[1613];
  //sp 4
  sp_rates[3] += fwd_rates[1613];
  //sp 100
  sp_rates[99] -= fwd_rates[1613];
  //sp 12
  sp_rates[11] += fwd_rates[1613];
  //sp 50
  sp_rates[49] -= fwd_rates[1613];

  //rxn 1614
  //sp 98
  sp_rates[97] += fwd_rates[1614];
  //sp 99
  sp_rates[98] -= fwd_rates[1614];
  //sp 4
  sp_rates[3] += fwd_rates[1614];
  //sp 12
  sp_rates[11] += fwd_rates[1614];
  //sp 50
  sp_rates[49] -= fwd_rates[1614];

  //rxn 1615
  //sp 4
  sp_rates[3] += fwd_rates[1615];
  //sp 12
  sp_rates[11] += fwd_rates[1615];
  //sp 50
  sp_rates[49] -= fwd_rates[1615];
  //sp 94
  sp_rates[93] += fwd_rates[1615];
  //sp 95
  sp_rates[94] -= fwd_rates[1615];

  //rxn 1616
  //sp 4
  sp_rates[3] += fwd_rates[1616];
  //sp 104
  sp_rates[103] += fwd_rates[1616];
  //sp 105
  sp_rates[104] -= fwd_rates[1616];
  //sp 12
  sp_rates[11] += fwd_rates[1616];
  //sp 50
  sp_rates[49] -= fwd_rates[1616];

  //rxn 1617
  //sp 4
  sp_rates[3] += fwd_rates[1617];
  //sp 107
  sp_rates[106] += fwd_rates[1617];
  //sp 12
  sp_rates[11] += fwd_rates[1617];
  //sp 109
  sp_rates[108] -= fwd_rates[1617];
  //sp 50
  sp_rates[49] -= fwd_rates[1617];

  //rxn 1618
  //sp 4
  sp_rates[3] += fwd_rates[1618];
  //sp 12
  sp_rates[11] += fwd_rates[1618];
  //sp 110
  sp_rates[109] += fwd_rates[1618];
  //sp 111
  sp_rates[110] -= fwd_rates[1618];
  //sp 50
  sp_rates[49] -= fwd_rates[1618];

  //rxn 1619
  //sp 4
  sp_rates[3] += fwd_rates[1619];
  //sp 6
  sp_rates[5] -= fwd_rates[1619];
  //sp 8
  sp_rates[7] += fwd_rates[1619];
  //sp 12
  sp_rates[11] += fwd_rates[1619];
  //sp 50
  sp_rates[49] -= fwd_rates[1619];

  //rxn 1620
  //sp 4
  sp_rates[3] += fwd_rates[1620];
  //sp 5
  sp_rates[4] -= fwd_rates[1620];
  //sp 10
  sp_rates[9] += fwd_rates[1620];
  //sp 12
  sp_rates[11] += fwd_rates[1620];
  //sp 50
  sp_rates[49] -= fwd_rates[1620];

  //rxn 1621
  //sp 4
  sp_rates[3] += fwd_rates[1621];
  //sp 7
  sp_rates[6] += fwd_rates[1621];
  //sp 8
  sp_rates[7] -= fwd_rates[1621];
  //sp 12
  sp_rates[11] += fwd_rates[1621];
  //sp 50
  sp_rates[49] -= fwd_rates[1621];

  //rxn 1622
  //sp 4
  sp_rates[3] += fwd_rates[1622];
  //sp 9
  sp_rates[8] += fwd_rates[1622];
  //sp 10
  sp_rates[9] -= fwd_rates[1622];
  //sp 12
  sp_rates[11] += fwd_rates[1622];
  //sp 50
  sp_rates[49] -= fwd_rates[1622];

  //rxn 1623
  //sp 4
  sp_rates[3] += fwd_rates[1623];
  //sp 12
  sp_rates[11] += fwd_rates[1623];
  //sp 50
  sp_rates[49] -= fwd_rates[1623];
  //sp 22
  sp_rates[21] += fwd_rates[1623];
  //sp 23
  sp_rates[22] -= fwd_rates[1623];

  //rxn 1624
  //sp 4
  sp_rates[3] += fwd_rates[1624];
  //sp 12
  sp_rates[11] += fwd_rates[1624];
  //sp 50
  sp_rates[49] -= fwd_rates[1624];
  //sp 22
  sp_rates[21] += fwd_rates[1624];
  //sp 24
  sp_rates[23] -= fwd_rates[1624];

  //rxn 1625
  //sp 4
  sp_rates[3] += fwd_rates[1625];
  //sp 12
  sp_rates[11] += fwd_rates[1625];
  //sp 50
  sp_rates[49] -= fwd_rates[1625];
  //sp 25
  sp_rates[24] += fwd_rates[1625];
  //sp 26
  sp_rates[25] -= fwd_rates[1625];

  //rxn 1626
  //sp 4
  sp_rates[3] += fwd_rates[1626];
  //sp 12
  sp_rates[11] += fwd_rates[1626];
  //sp 45
  sp_rates[44] += fwd_rates[1626];
  //sp 46
  sp_rates[45] -= fwd_rates[1626];
  //sp 50
  sp_rates[49] -= fwd_rates[1626];

  //rxn 1627
  //sp 4
  sp_rates[3] += fwd_rates[1627];
  //sp 12
  sp_rates[11] += fwd_rates[1627];
  //sp 45
  sp_rates[44] += fwd_rates[1627];
  //sp 47
  sp_rates[46] -= fwd_rates[1627];
  //sp 50
  sp_rates[49] -= fwd_rates[1627];

  //rxn 1628
  //sp 4
  sp_rates[3] += fwd_rates[1628];
  //sp 12
  sp_rates[11] += fwd_rates[1628];
  //sp 79
  sp_rates[78] += fwd_rates[1628];
  //sp 80
  sp_rates[79] -= fwd_rates[1628];
  //sp 50
  sp_rates[49] -= fwd_rates[1628];

  //rxn 1629
  //sp 4
  sp_rates[3] += fwd_rates[1629];
  //sp 12
  sp_rates[11] += fwd_rates[1629];
  //sp 48
  sp_rates[47] += fwd_rates[1629];
  //sp 49
  sp_rates[48] -= fwd_rates[1629];
  //sp 50
  sp_rates[49] -= fwd_rates[1629];

  //rxn 1630
  //sp 14
  sp_rates[13] += fwd_rates[1630];
  //sp 15
  sp_rates[14] -= fwd_rates[1630];
  //sp 57
  sp_rates[56] -= fwd_rates[1630];
  //sp 25
  sp_rates[24] += fwd_rates[1630];
  //sp 26
  sp_rates[25] += fwd_rates[1630];

  //rxn 1631
  //sp 34
  sp_rates[33] += fwd_rates[1631];
  //sp 35
  sp_rates[34] -= fwd_rates[1631];
  //sp 57
  sp_rates[56] -= fwd_rates[1631];
  //sp 25
  sp_rates[24] += fwd_rates[1631];
  //sp 26
  sp_rates[25] += fwd_rates[1631];

  //rxn 1632
  //sp 26
  sp_rates[25] += fwd_rates[1632];
  //sp 57
  sp_rates[56] -= fwd_rates[1632];
  //sp 25
  sp_rates[24] += fwd_rates[1632];
  //sp 30
  sp_rates[29] += fwd_rates[1632];
  //sp 31
  sp_rates[30] -= fwd_rates[1632];

  //rxn 1633
  //sp 57
  sp_rates[56] -= fwd_rates[1633];
  //sp 25
  sp_rates[24] += fwd_rates[1633];
  //sp 59
  sp_rates[58] += fwd_rates[1633];
  //sp 61
  sp_rates[60] -= fwd_rates[1633];
  //sp 26
  sp_rates[25] += fwd_rates[1633];

  //rxn 1634
  //sp 57
  sp_rates[56] -= fwd_rates[1634];
  //sp 25
  sp_rates[24] += fwd_rates[1634];
  //sp 59
  sp_rates[58] += fwd_rates[1634];
  //sp 60
  sp_rates[59] -= fwd_rates[1634];
  //sp 26
  sp_rates[25] += fwd_rates[1634];

  //rxn 1635
  //sp 57
  sp_rates[56] -= fwd_rates[1635];
  //sp 25
  sp_rates[24] += fwd_rates[1635];
  //sp 62
  sp_rates[61] += fwd_rates[1635];
  //sp 26
  sp_rates[25] += fwd_rates[1635];
  //sp 64
  sp_rates[63] -= fwd_rates[1635];

  //rxn 1636
  //sp 65
  sp_rates[64] -= fwd_rates[1636];
  //sp 57
  sp_rates[56] -= fwd_rates[1636];
  //sp 25
  sp_rates[24] += fwd_rates[1636];
  //sp 62
  sp_rates[61] += fwd_rates[1636];
  //sp 26
  sp_rates[25] += fwd_rates[1636];

  //rxn 1637
  //sp 26
  sp_rates[25] += fwd_rates[1637];
  //sp 57
  sp_rates[56] -= fwd_rates[1637];
  //sp 25
  sp_rates[24] += fwd_rates[1637];
  //sp 62
  sp_rates[61] += fwd_rates[1637];
  //sp 63
  sp_rates[62] -= fwd_rates[1637];

  //rxn 1638
  //sp 102
  sp_rates[101] += fwd_rates[1638];
  //sp 103
  sp_rates[102] -= fwd_rates[1638];
  //sp 57
  sp_rates[56] -= fwd_rates[1638];
  //sp 25
  sp_rates[24] += fwd_rates[1638];
  //sp 26
  sp_rates[25] += fwd_rates[1638];

  //rxn 1639
  //sp 97
  sp_rates[96] += fwd_rates[1639];
  //sp 100
  sp_rates[99] -= fwd_rates[1639];
  //sp 57
  sp_rates[56] -= fwd_rates[1639];
  //sp 25
  sp_rates[24] += fwd_rates[1639];
  //sp 26
  sp_rates[25] += fwd_rates[1639];

  //rxn 1640
  //sp 98
  sp_rates[97] += fwd_rates[1640];
  //sp 99
  sp_rates[98] -= fwd_rates[1640];
  //sp 57
  sp_rates[56] -= fwd_rates[1640];
  //sp 25
  sp_rates[24] += fwd_rates[1640];
  //sp 26
  sp_rates[25] += fwd_rates[1640];

  //rxn 1641
  //sp 26
  sp_rates[25] += fwd_rates[1641];
  //sp 57
  sp_rates[56] -= fwd_rates[1641];
  //sp 25
  sp_rates[24] += fwd_rates[1641];
  //sp 94
  sp_rates[93] += fwd_rates[1641];
  //sp 95
  sp_rates[94] -= fwd_rates[1641];

  //rxn 1642
  //sp 104
  sp_rates[103] += fwd_rates[1642];
  //sp 105
  sp_rates[104] -= fwd_rates[1642];
  //sp 57
  sp_rates[56] -= fwd_rates[1642];
  //sp 25
  sp_rates[24] += fwd_rates[1642];
  //sp 26
  sp_rates[25] += fwd_rates[1642];

  //rxn 1643
  //sp 107
  sp_rates[106] += fwd_rates[1643];
  //sp 109
  sp_rates[108] -= fwd_rates[1643];
  //sp 57
  sp_rates[56] -= fwd_rates[1643];
  //sp 25
  sp_rates[24] += fwd_rates[1643];
  //sp 26
  sp_rates[25] += fwd_rates[1643];

  //rxn 1644
  //sp 110
  sp_rates[109] += fwd_rates[1644];
  //sp 111
  sp_rates[110] -= fwd_rates[1644];
  //sp 57
  sp_rates[56] -= fwd_rates[1644];
  //sp 25
  sp_rates[24] += fwd_rates[1644];
  //sp 26
  sp_rates[25] += fwd_rates[1644];

  //rxn 1645
  //sp 6
  sp_rates[5] -= fwd_rates[1645];
  //sp 8
  sp_rates[7] += fwd_rates[1645];
  //sp 57
  sp_rates[56] -= fwd_rates[1645];
  //sp 25
  sp_rates[24] += fwd_rates[1645];
  //sp 26
  sp_rates[25] += fwd_rates[1645];

  //rxn 1646
  //sp 5
  sp_rates[4] -= fwd_rates[1646];
  //sp 10
  sp_rates[9] += fwd_rates[1646];
  //sp 57
  sp_rates[56] -= fwd_rates[1646];
  //sp 25
  sp_rates[24] += fwd_rates[1646];
  //sp 26
  sp_rates[25] += fwd_rates[1646];

  //rxn 1647
  //sp 7
  sp_rates[6] += fwd_rates[1647];
  //sp 8
  sp_rates[7] -= fwd_rates[1647];
  //sp 57
  sp_rates[56] -= fwd_rates[1647];
  //sp 25
  sp_rates[24] += fwd_rates[1647];
  //sp 26
  sp_rates[25] += fwd_rates[1647];

  //rxn 1648
  //sp 9
  sp_rates[8] += fwd_rates[1648];
  //sp 10
  sp_rates[9] -= fwd_rates[1648];
  //sp 57
  sp_rates[56] -= fwd_rates[1648];
  //sp 26
  sp_rates[25] += fwd_rates[1648];
  //sp 25
  sp_rates[24] += fwd_rates[1648];

  //rxn 1649
  //sp 22
  sp_rates[21] += fwd_rates[1649];
  //sp 23
  sp_rates[22] -= fwd_rates[1649];
  //sp 57
  sp_rates[56] -= fwd_rates[1649];
  //sp 25
  sp_rates[24] += fwd_rates[1649];
  //sp 26
  sp_rates[25] += fwd_rates[1649];

  //rxn 1650
  //sp 22
  sp_rates[21] += fwd_rates[1650];
  //sp 24
  sp_rates[23] -= fwd_rates[1650];
  //sp 57
  sp_rates[56] -= fwd_rates[1650];
  //sp 25
  sp_rates[24] += fwd_rates[1650];
  //sp 26
  sp_rates[25] += fwd_rates[1650];

  //rxn 1651
  //sp 57
  sp_rates[56] -= fwd_rates[1651];
  //sp 25
  sp_rates[24] += 2.0 * fwd_rates[1651];

  //rxn 1652
  //sp 45
  sp_rates[44] += fwd_rates[1652];
  //sp 46
  sp_rates[45] -= fwd_rates[1652];
  //sp 57
  sp_rates[56] -= fwd_rates[1652];
  //sp 25
  sp_rates[24] += fwd_rates[1652];
  //sp 26
  sp_rates[25] += fwd_rates[1652];

  //rxn 1653
  //sp 45
  sp_rates[44] += fwd_rates[1653];
  //sp 47
  sp_rates[46] -= fwd_rates[1653];
  //sp 57
  sp_rates[56] -= fwd_rates[1653];
  //sp 25
  sp_rates[24] += fwd_rates[1653];
  //sp 26
  sp_rates[25] += fwd_rates[1653];

  //rxn 1654
  //sp 79
  sp_rates[78] += fwd_rates[1654];
  //sp 80
  sp_rates[79] -= fwd_rates[1654];
  //sp 57
  sp_rates[56] -= fwd_rates[1654];
  //sp 25
  sp_rates[24] += fwd_rates[1654];
  //sp 26
  sp_rates[25] += fwd_rates[1654];

  //rxn 1655
  //sp 48
  sp_rates[47] += fwd_rates[1655];
  //sp 49
  sp_rates[48] -= fwd_rates[1655];
  //sp 25
  sp_rates[24] += fwd_rates[1655];
  //sp 26
  sp_rates[25] += fwd_rates[1655];
  //sp 57
  sp_rates[56] -= fwd_rates[1655];

  //rxn 1656
  //sp 57
  sp_rates[56] -= fwd_rates[1656];
  //sp 58
  sp_rates[57] += fwd_rates[1656];
  //sp 14
  sp_rates[13] += fwd_rates[1656];
  //sp 15
  sp_rates[14] -= fwd_rates[1656];

  //rxn 1657
  //sp 57
  sp_rates[56] -= fwd_rates[1657];
  //sp 34
  sp_rates[33] += fwd_rates[1657];
  //sp 35
  sp_rates[34] -= fwd_rates[1657];
  //sp 58
  sp_rates[57] += fwd_rates[1657];

  //rxn 1658
  //sp 57
  sp_rates[56] -= fwd_rates[1658];
  //sp 58
  sp_rates[57] += fwd_rates[1658];
  //sp 30
  sp_rates[29] += fwd_rates[1658];
  //sp 31
  sp_rates[30] -= fwd_rates[1658];

  //rxn 1659
  //sp 57
  sp_rates[56] -= fwd_rates[1659];
  //sp 58
  sp_rates[57] += fwd_rates[1659];
  //sp 59
  sp_rates[58] += fwd_rates[1659];
  //sp 61
  sp_rates[60] -= fwd_rates[1659];

  //rxn 1660
  //sp 57
  sp_rates[56] -= fwd_rates[1660];
  //sp 58
  sp_rates[57] += fwd_rates[1660];
  //sp 59
  sp_rates[58] += fwd_rates[1660];
  //sp 60
  sp_rates[59] -= fwd_rates[1660];

  //rxn 1661
  //sp 57
  sp_rates[56] -= fwd_rates[1661];
  //sp 58
  sp_rates[57] += fwd_rates[1661];
  //sp 62
  sp_rates[61] += fwd_rates[1661];
  //sp 64
  sp_rates[63] -= fwd_rates[1661];

  //rxn 1662
  //sp 57
  sp_rates[56] -= fwd_rates[1662];
  //sp 58
  sp_rates[57] += fwd_rates[1662];
  //sp 65
  sp_rates[64] -= fwd_rates[1662];
  //sp 62
  sp_rates[61] += fwd_rates[1662];

  //rxn 1663
  //sp 57
  sp_rates[56] -= fwd_rates[1663];
  //sp 58
  sp_rates[57] += fwd_rates[1663];
  //sp 62
  sp_rates[61] += fwd_rates[1663];
  //sp 63
  sp_rates[62] -= fwd_rates[1663];

  //rxn 1664
  //sp 57
  sp_rates[56] -= fwd_rates[1664];
  //sp 58
  sp_rates[57] += fwd_rates[1664];
  //sp 102
  sp_rates[101] += fwd_rates[1664];
  //sp 103
  sp_rates[102] -= fwd_rates[1664];

  //rxn 1665
  //sp 57
  sp_rates[56] -= fwd_rates[1665];
  //sp 58
  sp_rates[57] += fwd_rates[1665];
  //sp 100
  sp_rates[99] -= fwd_rates[1665];
  //sp 97
  sp_rates[96] += fwd_rates[1665];

  //rxn 1666
  //sp 57
  sp_rates[56] -= fwd_rates[1666];
  //sp 58
  sp_rates[57] += fwd_rates[1666];
  //sp 99
  sp_rates[98] -= fwd_rates[1666];
  //sp 98
  sp_rates[97] += fwd_rates[1666];

  //rxn 1667
  //sp 57
  sp_rates[56] -= fwd_rates[1667];
  //sp 58
  sp_rates[57] += fwd_rates[1667];
  //sp 94
  sp_rates[93] += fwd_rates[1667];
  //sp 95
  sp_rates[94] -= fwd_rates[1667];

  //rxn 1668
  //sp 57
  sp_rates[56] -= fwd_rates[1668];
  //sp 58
  sp_rates[57] += fwd_rates[1668];
  //sp 105
  sp_rates[104] -= fwd_rates[1668];
  //sp 104
  sp_rates[103] += fwd_rates[1668];

  //rxn 1669
  //sp 57
  sp_rates[56] -= fwd_rates[1669];
  //sp 58
  sp_rates[57] += fwd_rates[1669];
  //sp 107
  sp_rates[106] += fwd_rates[1669];
  //sp 109
  sp_rates[108] -= fwd_rates[1669];

  //rxn 1670
  //sp 57
  sp_rates[56] -= fwd_rates[1670];
  //sp 58
  sp_rates[57] += fwd_rates[1670];
  //sp 110
  sp_rates[109] += fwd_rates[1670];
  //sp 111
  sp_rates[110] -= fwd_rates[1670];

  //rxn 1671
  //sp 57
  sp_rates[56] -= fwd_rates[1671];
  //sp 58
  sp_rates[57] += fwd_rates[1671];
  //sp 6
  sp_rates[5] -= fwd_rates[1671];
  //sp 8
  sp_rates[7] += fwd_rates[1671];

  //rxn 1672
  //sp 57
  sp_rates[56] -= fwd_rates[1672];
  //sp 10
  sp_rates[9] += fwd_rates[1672];
  //sp 5
  sp_rates[4] -= fwd_rates[1672];
  //sp 58
  sp_rates[57] += fwd_rates[1672];

  //rxn 1673
  //sp 57
  sp_rates[56] -= fwd_rates[1673];
  //sp 58
  sp_rates[57] += fwd_rates[1673];
  //sp 7
  sp_rates[6] += fwd_rates[1673];
  //sp 8
  sp_rates[7] -= fwd_rates[1673];

  //rxn 1674
  //sp 57
  sp_rates[56] -= fwd_rates[1674];
  //sp 10
  sp_rates[9] -= fwd_rates[1674];
  //sp 9
  sp_rates[8] += fwd_rates[1674];
  //sp 58
  sp_rates[57] += fwd_rates[1674];

  //rxn 1675
  //sp 57
  sp_rates[56] -= fwd_rates[1675];
  //sp 58
  sp_rates[57] += fwd_rates[1675];
  //sp 22
  sp_rates[21] += fwd_rates[1675];
  //sp 23
  sp_rates[22] -= fwd_rates[1675];

  //rxn 1676
  //sp 57
  sp_rates[56] -= fwd_rates[1676];
  //sp 58
  sp_rates[57] += fwd_rates[1676];
  //sp 22
  sp_rates[21] += fwd_rates[1676];
  //sp 24
  sp_rates[23] -= fwd_rates[1676];

  //rxn 1677
  //sp 57
  sp_rates[56] -= fwd_rates[1677];
  //sp 26
  sp_rates[25] -= fwd_rates[1677];
  //sp 25
  sp_rates[24] += fwd_rates[1677];
  //sp 58
  sp_rates[57] += fwd_rates[1677];

  //rxn 1678
  //sp 57
  sp_rates[56] -= fwd_rates[1678];
  //sp 58
  sp_rates[57] += fwd_rates[1678];
  //sp 45
  sp_rates[44] += fwd_rates[1678];
  //sp 46
  sp_rates[45] -= fwd_rates[1678];

  //rxn 1679
  //sp 57
  sp_rates[56] -= fwd_rates[1679];
  //sp 58
  sp_rates[57] += fwd_rates[1679];
  //sp 45
  sp_rates[44] += fwd_rates[1679];
  //sp 47
  sp_rates[46] -= fwd_rates[1679];

  //rxn 1680
  //sp 57
  sp_rates[56] -= fwd_rates[1680];
  //sp 58
  sp_rates[57] += fwd_rates[1680];
  //sp 79
  sp_rates[78] += fwd_rates[1680];
  //sp 80
  sp_rates[79] -= fwd_rates[1680];

  //rxn 1681
  //sp 49
  sp_rates[48] -= fwd_rates[1681];
  //sp 58
  sp_rates[57] += fwd_rates[1681];
  //sp 57
  sp_rates[56] -= fwd_rates[1681];
  //sp 48
  sp_rates[47] += fwd_rates[1681];

  //rxn 1682
  //sp 20
  sp_rates[19] -= fwd_rates[1682];
  //sp 21
  sp_rates[20] += fwd_rates[1682];
  //sp 62
  sp_rates[61] += fwd_rates[1682];
  //sp 64
  sp_rates[63] -= fwd_rates[1682];

  //rxn 1683
  //sp 65
  sp_rates[64] -= fwd_rates[1683];
  //sp 20
  sp_rates[19] -= fwd_rates[1683];
  //sp 21
  sp_rates[20] += fwd_rates[1683];
  //sp 62
  sp_rates[61] += fwd_rates[1683];

  //rxn 1684
  //sp 20
  sp_rates[19] -= fwd_rates[1684];
  //sp 21
  sp_rates[20] += fwd_rates[1684];
  //sp 102
  sp_rates[101] += fwd_rates[1684];
  //sp 103
  sp_rates[102] -= fwd_rates[1684];

  //rxn 1685
  //sp 98
  sp_rates[97] += fwd_rates[1685];
  //sp 99
  sp_rates[98] -= fwd_rates[1685];
  //sp 20
  sp_rates[19] -= fwd_rates[1685];
  //sp 21
  sp_rates[20] += fwd_rates[1685];

  //rxn 1686
  //sp 105
  sp_rates[104] -= fwd_rates[1686];
  //sp 20
  sp_rates[19] -= fwd_rates[1686];
  //sp 21
  sp_rates[20] += fwd_rates[1686];
  //sp 104
  sp_rates[103] += fwd_rates[1686];

  //rxn 1687
  //sp 107
  sp_rates[106] += fwd_rates[1687];
  //sp 20
  sp_rates[19] -= fwd_rates[1687];
  //sp 109
  sp_rates[108] -= fwd_rates[1687];
  //sp 21
  sp_rates[20] += fwd_rates[1687];

  //rxn 1688
  //sp 20
  sp_rates[19] -= fwd_rates[1688];
  //sp 21
  sp_rates[20] += fwd_rates[1688];
  //sp 110
  sp_rates[109] += fwd_rates[1688];
  //sp 111
  sp_rates[110] -= fwd_rates[1688];

  //rxn 1689
  //sp 20
  sp_rates[19] -= fwd_rates[1689];
  //sp 21
  sp_rates[20] += fwd_rates[1689];
  //sp 6
  sp_rates[5] -= fwd_rates[1689];
  //sp 8
  sp_rates[7] += fwd_rates[1689];

  //rxn 1690
  //sp 20
  sp_rates[19] -= fwd_rates[1690];
  //sp 21
  sp_rates[20] += fwd_rates[1690];
  //sp 7
  sp_rates[6] += fwd_rates[1690];
  //sp 8
  sp_rates[7] -= fwd_rates[1690];

  //rxn 1691
  //sp 20
  sp_rates[19] -= fwd_rates[1691];
  //sp 21
  sp_rates[20] += fwd_rates[1691];
  //sp 22
  sp_rates[21] += fwd_rates[1691];
  //sp 23
  sp_rates[22] -= fwd_rates[1691];

  //rxn 1692
  //sp 20
  sp_rates[19] -= fwd_rates[1692];
  //sp 21
  sp_rates[20] += fwd_rates[1692];
  //sp 47
  sp_rates[46] -= fwd_rates[1692];
  //sp 45
  sp_rates[44] += fwd_rates[1692];

  //rxn 1693
  //sp 49
  sp_rates[48] -= fwd_rates[1693];
  //sp 20
  sp_rates[19] -= fwd_rates[1693];
  //sp 21
  sp_rates[20] += fwd_rates[1693];
  //sp 48
  sp_rates[47] += fwd_rates[1693];

  //rxn 1694
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1694];
  //sp 3
  sp_rates[2] += fwd_rates[1694];
  //sp 4
  sp_rates[3] -= fwd_rates[1694];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1694];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1694];
  //sp 32
  sp_rates[31] -= fwd_rates[1694];

  //rxn 1695
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1695];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1695];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1695];
  //sp 14
  sp_rates[13] += fwd_rates[1695];
  //sp 15
  sp_rates[14] -= fwd_rates[1695];
  //sp 32
  sp_rates[31] -= fwd_rates[1695];

  //rxn 1696
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1696];
  //sp 34
  sp_rates[33] += fwd_rates[1696];
  //sp 35
  sp_rates[34] -= fwd_rates[1696];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1696];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1696];
  //sp 32
  sp_rates[31] -= fwd_rates[1696];

  //rxn 1697
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1697];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1697];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1697];
  //sp 30
  sp_rates[29] += fwd_rates[1697];
  //sp 31
  sp_rates[30] -= fwd_rates[1697];
  //sp 32
  sp_rates[31] -= fwd_rates[1697];

  //rxn 1698
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1698];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1698];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1698];
  //sp 59
  sp_rates[58] += fwd_rates[1698];
  //sp 61
  sp_rates[60] -= fwd_rates[1698];
  //sp 32
  sp_rates[31] -= fwd_rates[1698];

  //rxn 1699
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1699];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1699];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1699];
  //sp 59
  sp_rates[58] += fwd_rates[1699];
  //sp 60
  sp_rates[59] -= fwd_rates[1699];
  //sp 32
  sp_rates[31] -= fwd_rates[1699];

  //rxn 1700
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1700];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1700];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1700];
  //sp 32
  sp_rates[31] -= fwd_rates[1700];
  //sp 62
  sp_rates[61] += fwd_rates[1700];
  //sp 64
  sp_rates[63] -= fwd_rates[1700];

  //rxn 1701
  //sp 65
  sp_rates[64] -= fwd_rates[1701];
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1701];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1701];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1701];
  //sp 62
  sp_rates[61] += fwd_rates[1701];
  //sp 32
  sp_rates[31] -= fwd_rates[1701];

  //rxn 1702
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1702];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1702];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1702];
  //sp 62
  sp_rates[61] += fwd_rates[1702];
  //sp 63
  sp_rates[62] -= fwd_rates[1702];
  //sp 32
  sp_rates[31] -= fwd_rates[1702];

  //rxn 1703
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1703];
  //sp 102
  sp_rates[101] += fwd_rates[1703];
  //sp 103
  sp_rates[102] -= fwd_rates[1703];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1703];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1703];
  //sp 32
  sp_rates[31] -= fwd_rates[1703];

  //rxn 1704
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1704];
  //sp 97
  sp_rates[96] += fwd_rates[1704];
  //sp 100
  sp_rates[99] -= fwd_rates[1704];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1704];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1704];
  //sp 32
  sp_rates[31] -= fwd_rates[1704];

  //rxn 1705
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1705];
  //sp 98
  sp_rates[97] += fwd_rates[1705];
  //sp 99
  sp_rates[98] -= fwd_rates[1705];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1705];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1705];
  //sp 32
  sp_rates[31] -= fwd_rates[1705];

  //rxn 1706
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1706];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1706];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1706];
  //sp 94
  sp_rates[93] += fwd_rates[1706];
  //sp 95
  sp_rates[94] -= fwd_rates[1706];
  //sp 32
  sp_rates[31] -= fwd_rates[1706];

  //rxn 1707
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1707];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1707];
  //sp 105
  sp_rates[104] -= fwd_rates[1707];
  //sp 104
  sp_rates[103] += fwd_rates[1707];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1707];
  //sp 32
  sp_rates[31] -= fwd_rates[1707];

  //rxn 1708
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1708];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1708];
  //sp 107
  sp_rates[106] += fwd_rates[1708];
  //sp 109
  sp_rates[108] -= fwd_rates[1708];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1708];
  //sp 32
  sp_rates[31] -= fwd_rates[1708];

  //rxn 1709
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1709];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1709];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1709];
  //sp 110
  sp_rates[109] += fwd_rates[1709];
  //sp 111
  sp_rates[110] -= fwd_rates[1709];
  //sp 32
  sp_rates[31] -= fwd_rates[1709];

  //rxn 1710
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1710];
  //sp 6
  sp_rates[5] -= fwd_rates[1710];
  //sp 8
  sp_rates[7] += 1.5 * fwd_rates[1710];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1710];
  //sp 32
  sp_rates[31] -= fwd_rates[1710];

  //rxn 1711
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1711];
  //sp 5
  sp_rates[4] -= fwd_rates[1711];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1711];
  //sp 10
  sp_rates[9] += fwd_rates[1711];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1711];
  //sp 32
  sp_rates[31] -= fwd_rates[1711];

  //rxn 1712
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1712];
  //sp 7
  sp_rates[6] += fwd_rates[1712];
  //sp 8
  sp_rates[7] -= 0.5 * fwd_rates[1712];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1712];
  //sp 32
  sp_rates[31] -= fwd_rates[1712];

  //rxn 1713
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1713];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1713];
  //sp 9
  sp_rates[8] += fwd_rates[1713];
  //sp 10
  sp_rates[9] -= fwd_rates[1713];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1713];
  //sp 32
  sp_rates[31] -= fwd_rates[1713];

  //rxn 1714
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1714];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1714];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1714];
  //sp 22
  sp_rates[21] += fwd_rates[1714];
  //sp 23
  sp_rates[22] -= fwd_rates[1714];
  //sp 32
  sp_rates[31] -= fwd_rates[1714];

  //rxn 1715
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1715];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1715];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1715];
  //sp 22
  sp_rates[21] += fwd_rates[1715];
  //sp 24
  sp_rates[23] -= fwd_rates[1715];
  //sp 32
  sp_rates[31] -= fwd_rates[1715];

  //rxn 1716
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1716];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1716];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1716];
  //sp 25
  sp_rates[24] += fwd_rates[1716];
  //sp 26
  sp_rates[25] -= fwd_rates[1716];
  //sp 32
  sp_rates[31] -= fwd_rates[1716];

  //rxn 1717
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1717];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1717];
  //sp 45
  sp_rates[44] += 1.5 * fwd_rates[1717];
  //sp 46
  sp_rates[45] -= fwd_rates[1717];
  //sp 32
  sp_rates[31] -= fwd_rates[1717];

  //rxn 1718
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1718];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1718];
  //sp 45
  sp_rates[44] += 1.5 * fwd_rates[1718];
  //sp 47
  sp_rates[46] -= fwd_rates[1718];
  //sp 32
  sp_rates[31] -= fwd_rates[1718];

  //rxn 1719
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1719];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1719];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1719];
  //sp 79
  sp_rates[78] += fwd_rates[1719];
  //sp 80
  sp_rates[79] -= fwd_rates[1719];
  //sp 32
  sp_rates[31] -= fwd_rates[1719];

  //rxn 1720
  //sp 33
  sp_rates[32] += 0.5 * fwd_rates[1720];
  //sp 8
  sp_rates[7] += 0.5 * fwd_rates[1720];
  //sp 45
  sp_rates[44] += 0.5 * fwd_rates[1720];
  //sp 48
  sp_rates[47] += fwd_rates[1720];
  //sp 49
  sp_rates[48] -= fwd_rates[1720];
  //sp 32
  sp_rates[31] -= fwd_rates[1720];

  //rxn 1721
  //sp 3
  sp_rates[2] += fwd_rates[1721];
  //sp 4
  sp_rates[3] -= fwd_rates[1721];
  //sp 35
  sp_rates[34] += fwd_rates[1721];
  //sp 11
  sp_rates[10] += fwd_rates[1721];
  //sp 85
  sp_rates[84] -= fwd_rates[1721];

  //rxn 1722
  //sp 35
  sp_rates[34] += fwd_rates[1722];
  //sp 11
  sp_rates[10] += fwd_rates[1722];
  //sp 14
  sp_rates[13] += fwd_rates[1722];
  //sp 15
  sp_rates[14] -= fwd_rates[1722];
  //sp 85
  sp_rates[84] -= fwd_rates[1722];

  //rxn 1723
  //sp 34
  sp_rates[33] += fwd_rates[1723];
  //sp 11
  sp_rates[10] += fwd_rates[1723];
  //sp 85
  sp_rates[84] -= fwd_rates[1723];

  //rxn 1724
  //sp 35
  sp_rates[34] += fwd_rates[1724];
  //sp 11
  sp_rates[10] += fwd_rates[1724];
  //sp 85
  sp_rates[84] -= fwd_rates[1724];
  //sp 30
  sp_rates[29] += fwd_rates[1724];
  //sp 31
  sp_rates[30] -= fwd_rates[1724];

  //rxn 1725
  //sp 35
  sp_rates[34] += fwd_rates[1725];
  //sp 11
  sp_rates[10] += fwd_rates[1725];
  //sp 85
  sp_rates[84] -= fwd_rates[1725];
  //sp 59
  sp_rates[58] += fwd_rates[1725];
  //sp 61
  sp_rates[60] -= fwd_rates[1725];

  //rxn 1726
  //sp 35
  sp_rates[34] += fwd_rates[1726];
  //sp 11
  sp_rates[10] += fwd_rates[1726];
  //sp 85
  sp_rates[84] -= fwd_rates[1726];
  //sp 59
  sp_rates[58] += fwd_rates[1726];
  //sp 60
  sp_rates[59] -= fwd_rates[1726];

  //rxn 1727
  //sp 35
  sp_rates[34] += fwd_rates[1727];
  //sp 11
  sp_rates[10] += fwd_rates[1727];
  //sp 85
  sp_rates[84] -= fwd_rates[1727];
  //sp 62
  sp_rates[61] += fwd_rates[1727];
  //sp 64
  sp_rates[63] -= fwd_rates[1727];

  //rxn 1728
  //sp 65
  sp_rates[64] -= fwd_rates[1728];
  //sp 35
  sp_rates[34] += fwd_rates[1728];
  //sp 11
  sp_rates[10] += fwd_rates[1728];
  //sp 85
  sp_rates[84] -= fwd_rates[1728];
  //sp 62
  sp_rates[61] += fwd_rates[1728];

  //rxn 1729
  //sp 35
  sp_rates[34] += fwd_rates[1729];
  //sp 11
  sp_rates[10] += fwd_rates[1729];
  //sp 85
  sp_rates[84] -= fwd_rates[1729];
  //sp 62
  sp_rates[61] += fwd_rates[1729];
  //sp 63
  sp_rates[62] -= fwd_rates[1729];

  //rxn 1730
  //sp 35
  sp_rates[34] += fwd_rates[1730];
  //sp 102
  sp_rates[101] += fwd_rates[1730];
  //sp 103
  sp_rates[102] -= fwd_rates[1730];
  //sp 11
  sp_rates[10] += fwd_rates[1730];
  //sp 85
  sp_rates[84] -= fwd_rates[1730];

  //rxn 1731
  //sp 97
  sp_rates[96] += fwd_rates[1731];
  //sp 35
  sp_rates[34] += fwd_rates[1731];
  //sp 100
  sp_rates[99] -= fwd_rates[1731];
  //sp 11
  sp_rates[10] += fwd_rates[1731];
  //sp 85
  sp_rates[84] -= fwd_rates[1731];

  //rxn 1732
  //sp 98
  sp_rates[97] += fwd_rates[1732];
  //sp 99
  sp_rates[98] -= fwd_rates[1732];
  //sp 35
  sp_rates[34] += fwd_rates[1732];
  //sp 11
  sp_rates[10] += fwd_rates[1732];
  //sp 85
  sp_rates[84] -= fwd_rates[1732];

  //rxn 1733
  //sp 35
  sp_rates[34] += fwd_rates[1733];
  //sp 11
  sp_rates[10] += fwd_rates[1733];
  //sp 85
  sp_rates[84] -= fwd_rates[1733];
  //sp 94
  sp_rates[93] += fwd_rates[1733];
  //sp 95
  sp_rates[94] -= fwd_rates[1733];

  //rxn 1734
  //sp 35
  sp_rates[34] += fwd_rates[1734];
  //sp 104
  sp_rates[103] += fwd_rates[1734];
  //sp 105
  sp_rates[104] -= fwd_rates[1734];
  //sp 11
  sp_rates[10] += fwd_rates[1734];
  //sp 85
  sp_rates[84] -= fwd_rates[1734];

  //rxn 1735
  //sp 35
  sp_rates[34] += fwd_rates[1735];
  //sp 11
  sp_rates[10] += fwd_rates[1735];
  //sp 107
  sp_rates[106] += fwd_rates[1735];
  //sp 109
  sp_rates[108] -= fwd_rates[1735];
  //sp 85
  sp_rates[84] -= fwd_rates[1735];

  //rxn 1736
  //sp 35
  sp_rates[34] += fwd_rates[1736];
  //sp 11
  sp_rates[10] += fwd_rates[1736];
  //sp 110
  sp_rates[109] += fwd_rates[1736];
  //sp 111
  sp_rates[110] -= fwd_rates[1736];
  //sp 85
  sp_rates[84] -= fwd_rates[1736];

  //rxn 1737
  //sp 35
  sp_rates[34] += fwd_rates[1737];
  //sp 6
  sp_rates[5] -= fwd_rates[1737];
  //sp 8
  sp_rates[7] += fwd_rates[1737];
  //sp 11
  sp_rates[10] += fwd_rates[1737];
  //sp 85
  sp_rates[84] -= fwd_rates[1737];

  //rxn 1738
  //sp 35
  sp_rates[34] += fwd_rates[1738];
  //sp 5
  sp_rates[4] -= fwd_rates[1738];
  //sp 10
  sp_rates[9] += fwd_rates[1738];
  //sp 11
  sp_rates[10] += fwd_rates[1738];
  //sp 85
  sp_rates[84] -= fwd_rates[1738];

  //rxn 1739
  //sp 35
  sp_rates[34] += fwd_rates[1739];
  //sp 9
  sp_rates[8] += fwd_rates[1739];
  //sp 10
  sp_rates[9] -= fwd_rates[1739];
  //sp 11
  sp_rates[10] += fwd_rates[1739];
  //sp 85
  sp_rates[84] -= fwd_rates[1739];

  //rxn 1740
  //sp 35
  sp_rates[34] += fwd_rates[1740];
  //sp 11
  sp_rates[10] += fwd_rates[1740];
  //sp 85
  sp_rates[84] -= fwd_rates[1740];
  //sp 22
  sp_rates[21] += fwd_rates[1740];
  //sp 23
  sp_rates[22] -= fwd_rates[1740];

  //rxn 1741
  //sp 35
  sp_rates[34] += fwd_rates[1741];
  //sp 11
  sp_rates[10] += fwd_rates[1741];
  //sp 85
  sp_rates[84] -= fwd_rates[1741];
  //sp 22
  sp_rates[21] += fwd_rates[1741];
  //sp 24
  sp_rates[23] -= fwd_rates[1741];

  //rxn 1742
  //sp 35
  sp_rates[34] += fwd_rates[1742];
  //sp 11
  sp_rates[10] += fwd_rates[1742];
  //sp 85
  sp_rates[84] -= fwd_rates[1742];
  //sp 25
  sp_rates[24] += fwd_rates[1742];
  //sp 26
  sp_rates[25] -= fwd_rates[1742];

  //rxn 1743
  //sp 35
  sp_rates[34] += fwd_rates[1743];
  //sp 11
  sp_rates[10] += fwd_rates[1743];
  //sp 45
  sp_rates[44] += fwd_rates[1743];
  //sp 46
  sp_rates[45] -= fwd_rates[1743];
  //sp 85
  sp_rates[84] -= fwd_rates[1743];

  //rxn 1744
  //sp 35
  sp_rates[34] += fwd_rates[1744];
  //sp 11
  sp_rates[10] += fwd_rates[1744];
  //sp 45
  sp_rates[44] += fwd_rates[1744];
  //sp 47
  sp_rates[46] -= fwd_rates[1744];
  //sp 85
  sp_rates[84] -= fwd_rates[1744];

  //rxn 1745
  //sp 35
  sp_rates[34] += fwd_rates[1745];
  //sp 11
  sp_rates[10] += fwd_rates[1745];
  //sp 79
  sp_rates[78] += fwd_rates[1745];
  //sp 80
  sp_rates[79] -= fwd_rates[1745];
  //sp 85
  sp_rates[84] -= fwd_rates[1745];

  //rxn 1746
  //sp 35
  sp_rates[34] += fwd_rates[1746];
  //sp 11
  sp_rates[10] += fwd_rates[1746];
  //sp 48
  sp_rates[47] += fwd_rates[1746];
  //sp 49
  sp_rates[48] -= fwd_rates[1746];
  //sp 85
  sp_rates[84] -= fwd_rates[1746];

  //rxn 1747
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1747];
  //sp 3
  sp_rates[2] += fwd_rates[1747];
  //sp 67
  sp_rates[66] -= fwd_rates[1747];
  //sp 4
  sp_rates[3] -= fwd_rates[1747];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1747];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1747];

  //rxn 1748
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1748];
  //sp 67
  sp_rates[66] -= fwd_rates[1748];
  //sp 14
  sp_rates[13] += fwd_rates[1748];
  //sp 15
  sp_rates[14] -= fwd_rates[1748];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1748];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1748];

  //rxn 1749
  //sp 34
  sp_rates[33] += 1.3 * fwd_rates[1749];
  //sp 35
  sp_rates[34] -= fwd_rates[1749];
  //sp 67
  sp_rates[66] -= fwd_rates[1749];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1749];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1749];

  //rxn 1750
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1750];
  //sp 67
  sp_rates[66] -= fwd_rates[1750];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1750];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1750];
  //sp 30
  sp_rates[29] += fwd_rates[1750];
  //sp 31
  sp_rates[30] -= fwd_rates[1750];

  //rxn 1751
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1751];
  //sp 67
  sp_rates[66] -= fwd_rates[1751];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1751];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1751];
  //sp 59
  sp_rates[58] += fwd_rates[1751];
  //sp 61
  sp_rates[60] -= fwd_rates[1751];

  //rxn 1752
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1752];
  //sp 67
  sp_rates[66] -= fwd_rates[1752];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1752];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1752];
  //sp 59
  sp_rates[58] += fwd_rates[1752];
  //sp 60
  sp_rates[59] -= fwd_rates[1752];

  //rxn 1753
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1753];
  //sp 67
  sp_rates[66] -= fwd_rates[1753];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1753];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1753];
  //sp 62
  sp_rates[61] += fwd_rates[1753];
  //sp 64
  sp_rates[63] -= fwd_rates[1753];

  //rxn 1754
  //sp 65
  sp_rates[64] -= fwd_rates[1754];
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1754];
  //sp 67
  sp_rates[66] -= fwd_rates[1754];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1754];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1754];
  //sp 62
  sp_rates[61] += fwd_rates[1754];

  //rxn 1755
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1755];
  //sp 67
  sp_rates[66] -= fwd_rates[1755];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1755];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1755];
  //sp 62
  sp_rates[61] += fwd_rates[1755];
  //sp 63
  sp_rates[62] -= fwd_rates[1755];

  //rxn 1756
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1756];
  //sp 67
  sp_rates[66] -= fwd_rates[1756];
  //sp 102
  sp_rates[101] += fwd_rates[1756];
  //sp 103
  sp_rates[102] -= fwd_rates[1756];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1756];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1756];

  //rxn 1757
  //sp 97
  sp_rates[96] += fwd_rates[1757];
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1757];
  //sp 67
  sp_rates[66] -= fwd_rates[1757];
  //sp 100
  sp_rates[99] -= fwd_rates[1757];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1757];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1757];

  //rxn 1758
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1758];
  //sp 67
  sp_rates[66] -= fwd_rates[1758];
  //sp 99
  sp_rates[98] -= fwd_rates[1758];
  //sp 98
  sp_rates[97] += fwd_rates[1758];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1758];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1758];

  //rxn 1759
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1759];
  //sp 67
  sp_rates[66] -= fwd_rates[1759];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1759];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1759];
  //sp 94
  sp_rates[93] += fwd_rates[1759];
  //sp 95
  sp_rates[94] -= fwd_rates[1759];

  //rxn 1760
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1760];
  //sp 67
  sp_rates[66] -= fwd_rates[1760];
  //sp 104
  sp_rates[103] += fwd_rates[1760];
  //sp 105
  sp_rates[104] -= fwd_rates[1760];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1760];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1760];

  //rxn 1761
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1761];
  //sp 67
  sp_rates[66] -= fwd_rates[1761];
  //sp 107
  sp_rates[106] += fwd_rates[1761];
  //sp 109
  sp_rates[108] -= fwd_rates[1761];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1761];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1761];

  //rxn 1762
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1762];
  //sp 67
  sp_rates[66] -= fwd_rates[1762];
  //sp 110
  sp_rates[109] += fwd_rates[1762];
  //sp 111
  sp_rates[110] -= fwd_rates[1762];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1762];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1762];

  //rxn 1763
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1763];
  //sp 67
  sp_rates[66] -= fwd_rates[1763];
  //sp 6
  sp_rates[5] -= fwd_rates[1763];
  //sp 8
  sp_rates[7] += fwd_rates[1763];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1763];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1763];

  //rxn 1764
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1764];
  //sp 67
  sp_rates[66] -= fwd_rates[1764];
  //sp 5
  sp_rates[4] -= fwd_rates[1764];
  //sp 10
  sp_rates[9] += fwd_rates[1764];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1764];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1764];

  //rxn 1765
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1765];
  //sp 67
  sp_rates[66] -= fwd_rates[1765];
  //sp 7
  sp_rates[6] += fwd_rates[1765];
  //sp 8
  sp_rates[7] -= fwd_rates[1765];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1765];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1765];

  //rxn 1766
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1766];
  //sp 67
  sp_rates[66] -= fwd_rates[1766];
  //sp 9
  sp_rates[8] += fwd_rates[1766];
  //sp 10
  sp_rates[9] -= fwd_rates[1766];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1766];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1766];

  //rxn 1767
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1767];
  //sp 67
  sp_rates[66] -= fwd_rates[1767];
  //sp 22
  sp_rates[21] += fwd_rates[1767];
  //sp 23
  sp_rates[22] -= fwd_rates[1767];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1767];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1767];

  //rxn 1768
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1768];
  //sp 67
  sp_rates[66] -= fwd_rates[1768];
  //sp 22
  sp_rates[21] += fwd_rates[1768];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1768];
  //sp 24
  sp_rates[23] -= fwd_rates[1768];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1768];

  //rxn 1769
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1769];
  //sp 67
  sp_rates[66] -= fwd_rates[1769];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1769];
  //sp 25
  sp_rates[24] += fwd_rates[1769];
  //sp 26
  sp_rates[25] -= 0.7 * fwd_rates[1769];

  //rxn 1770
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1770];
  //sp 67
  sp_rates[66] -= fwd_rates[1770];
  //sp 45
  sp_rates[44] += fwd_rates[1770];
  //sp 46
  sp_rates[45] -= fwd_rates[1770];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1770];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1770];

  //rxn 1771
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1771];
  //sp 67
  sp_rates[66] -= fwd_rates[1771];
  //sp 45
  sp_rates[44] += fwd_rates[1771];
  //sp 47
  sp_rates[46] -= fwd_rates[1771];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1771];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1771];

  //rxn 1772
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1772];
  //sp 67
  sp_rates[66] -= fwd_rates[1772];
  //sp 79
  sp_rates[78] += fwd_rates[1772];
  //sp 80
  sp_rates[79] -= fwd_rates[1772];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1772];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1772];

  //rxn 1773
  //sp 34
  sp_rates[33] += 0.3 * fwd_rates[1773];
  //sp 67
  sp_rates[66] -= fwd_rates[1773];
  //sp 48
  sp_rates[47] += fwd_rates[1773];
  //sp 49
  sp_rates[48] -= fwd_rates[1773];
  //sp 87
  sp_rates[86] += 0.7 * fwd_rates[1773];
  //sp 26
  sp_rates[25] += 0.3 * fwd_rates[1773];

  //rxn 1774
  //sp 3
  sp_rates[2] += fwd_rates[1774];
  //sp 67
  sp_rates[66] -= fwd_rates[1774];
  //sp 85
  sp_rates[84] += fwd_rates[1774];

  //rxn 1775
  //sp 67
  sp_rates[66] -= fwd_rates[1775];
  //sp 4
  sp_rates[3] += fwd_rates[1775];
  //sp 14
  sp_rates[13] += fwd_rates[1775];
  //sp 15
  sp_rates[14] -= fwd_rates[1775];
  //sp 85
  sp_rates[84] += fwd_rates[1775];

  //rxn 1776
  //sp 34
  sp_rates[33] += fwd_rates[1776];
  //sp 35
  sp_rates[34] -= fwd_rates[1776];
  //sp 4
  sp_rates[3] += fwd_rates[1776];
  //sp 67
  sp_rates[66] -= fwd_rates[1776];
  //sp 85
  sp_rates[84] += fwd_rates[1776];

  //rxn 1777
  //sp 67
  sp_rates[66] -= fwd_rates[1777];
  //sp 4
  sp_rates[3] += fwd_rates[1777];
  //sp 85
  sp_rates[84] += fwd_rates[1777];
  //sp 30
  sp_rates[29] += fwd_rates[1777];
  //sp 31
  sp_rates[30] -= fwd_rates[1777];

  //rxn 1778
  //sp 67
  sp_rates[66] -= fwd_rates[1778];
  //sp 4
  sp_rates[3] += fwd_rates[1778];
  //sp 85
  sp_rates[84] += fwd_rates[1778];
  //sp 59
  sp_rates[58] += fwd_rates[1778];
  //sp 61
  sp_rates[60] -= fwd_rates[1778];

  //rxn 1779
  //sp 67
  sp_rates[66] -= fwd_rates[1779];
  //sp 4
  sp_rates[3] += fwd_rates[1779];
  //sp 85
  sp_rates[84] += fwd_rates[1779];
  //sp 59
  sp_rates[58] += fwd_rates[1779];
  //sp 60
  sp_rates[59] -= fwd_rates[1779];

  //rxn 1780
  //sp 67
  sp_rates[66] -= fwd_rates[1780];
  //sp 4
  sp_rates[3] += fwd_rates[1780];
  //sp 85
  sp_rates[84] += fwd_rates[1780];
  //sp 62
  sp_rates[61] += fwd_rates[1780];
  //sp 64
  sp_rates[63] -= fwd_rates[1780];

  //rxn 1781
  //sp 65
  sp_rates[64] -= fwd_rates[1781];
  //sp 67
  sp_rates[66] -= fwd_rates[1781];
  //sp 4
  sp_rates[3] += fwd_rates[1781];
  //sp 85
  sp_rates[84] += fwd_rates[1781];
  //sp 62
  sp_rates[61] += fwd_rates[1781];

  //rxn 1782
  //sp 67
  sp_rates[66] -= fwd_rates[1782];
  //sp 4
  sp_rates[3] += fwd_rates[1782];
  //sp 85
  sp_rates[84] += fwd_rates[1782];
  //sp 62
  sp_rates[61] += fwd_rates[1782];
  //sp 63
  sp_rates[62] -= fwd_rates[1782];

  //rxn 1783
  //sp 67
  sp_rates[66] -= fwd_rates[1783];
  //sp 4
  sp_rates[3] += fwd_rates[1783];
  //sp 102
  sp_rates[101] += fwd_rates[1783];
  //sp 103
  sp_rates[102] -= fwd_rates[1783];
  //sp 85
  sp_rates[84] += fwd_rates[1783];

  //rxn 1784
  //sp 97
  sp_rates[96] += fwd_rates[1784];
  //sp 67
  sp_rates[66] -= fwd_rates[1784];
  //sp 4
  sp_rates[3] += fwd_rates[1784];
  //sp 100
  sp_rates[99] -= fwd_rates[1784];
  //sp 85
  sp_rates[84] += fwd_rates[1784];

  //rxn 1785
  //sp 98
  sp_rates[97] += fwd_rates[1785];
  //sp 67
  sp_rates[66] -= fwd_rates[1785];
  //sp 4
  sp_rates[3] += fwd_rates[1785];
  //sp 99
  sp_rates[98] -= fwd_rates[1785];
  //sp 85
  sp_rates[84] += fwd_rates[1785];

  //rxn 1786
  //sp 67
  sp_rates[66] -= fwd_rates[1786];
  //sp 4
  sp_rates[3] += fwd_rates[1786];
  //sp 85
  sp_rates[84] += fwd_rates[1786];
  //sp 94
  sp_rates[93] += fwd_rates[1786];
  //sp 95
  sp_rates[94] -= fwd_rates[1786];

  //rxn 1787
  //sp 67
  sp_rates[66] -= fwd_rates[1787];
  //sp 4
  sp_rates[3] += fwd_rates[1787];
  //sp 104
  sp_rates[103] += fwd_rates[1787];
  //sp 105
  sp_rates[104] -= fwd_rates[1787];
  //sp 85
  sp_rates[84] += fwd_rates[1787];

  //rxn 1788
  //sp 67
  sp_rates[66] -= fwd_rates[1788];
  //sp 4
  sp_rates[3] += fwd_rates[1788];
  //sp 107
  sp_rates[106] += fwd_rates[1788];
  //sp 109
  sp_rates[108] -= fwd_rates[1788];
  //sp 85
  sp_rates[84] += fwd_rates[1788];

  //rxn 1789
  //sp 67
  sp_rates[66] -= fwd_rates[1789];
  //sp 4
  sp_rates[3] += fwd_rates[1789];
  //sp 110
  sp_rates[109] += fwd_rates[1789];
  //sp 111
  sp_rates[110] -= fwd_rates[1789];
  //sp 85
  sp_rates[84] += fwd_rates[1789];

  //rxn 1790
  //sp 67
  sp_rates[66] -= fwd_rates[1790];
  //sp 4
  sp_rates[3] += fwd_rates[1790];
  //sp 6
  sp_rates[5] -= fwd_rates[1790];
  //sp 8
  sp_rates[7] += fwd_rates[1790];
  //sp 85
  sp_rates[84] += fwd_rates[1790];

  //rxn 1791
  //sp 67
  sp_rates[66] -= fwd_rates[1791];
  //sp 4
  sp_rates[3] += fwd_rates[1791];
  //sp 5
  sp_rates[4] -= fwd_rates[1791];
  //sp 10
  sp_rates[9] += fwd_rates[1791];
  //sp 85
  sp_rates[84] += fwd_rates[1791];

  //rxn 1792
  //sp 67
  sp_rates[66] -= fwd_rates[1792];
  //sp 4
  sp_rates[3] += fwd_rates[1792];
  //sp 7
  sp_rates[6] += fwd_rates[1792];
  //sp 8
  sp_rates[7] -= fwd_rates[1792];
  //sp 85
  sp_rates[84] += fwd_rates[1792];

  //rxn 1793
  //sp 67
  sp_rates[66] -= fwd_rates[1793];
  //sp 4
  sp_rates[3] += fwd_rates[1793];
  //sp 9
  sp_rates[8] += fwd_rates[1793];
  //sp 10
  sp_rates[9] -= fwd_rates[1793];
  //sp 85
  sp_rates[84] += fwd_rates[1793];

  //rxn 1794
  //sp 67
  sp_rates[66] -= fwd_rates[1794];
  //sp 4
  sp_rates[3] += fwd_rates[1794];
  //sp 85
  sp_rates[84] += fwd_rates[1794];
  //sp 22
  sp_rates[21] += fwd_rates[1794];
  //sp 23
  sp_rates[22] -= fwd_rates[1794];

  //rxn 1795
  //sp 67
  sp_rates[66] -= fwd_rates[1795];
  //sp 4
  sp_rates[3] += fwd_rates[1795];
  //sp 85
  sp_rates[84] += fwd_rates[1795];
  //sp 22
  sp_rates[21] += fwd_rates[1795];
  //sp 24
  sp_rates[23] -= fwd_rates[1795];

  //rxn 1796
  //sp 67
  sp_rates[66] -= fwd_rates[1796];
  //sp 4
  sp_rates[3] += fwd_rates[1796];
  //sp 85
  sp_rates[84] += fwd_rates[1796];
  //sp 25
  sp_rates[24] += fwd_rates[1796];
  //sp 26
  sp_rates[25] -= fwd_rates[1796];

  //rxn 1797
  //sp 67
  sp_rates[66] -= fwd_rates[1797];
  //sp 4
  sp_rates[3] += fwd_rates[1797];
  //sp 45
  sp_rates[44] += fwd_rates[1797];
  //sp 46
  sp_rates[45] -= fwd_rates[1797];
  //sp 85
  sp_rates[84] += fwd_rates[1797];

  //rxn 1798
  //sp 67
  sp_rates[66] -= fwd_rates[1798];
  //sp 4
  sp_rates[3] += fwd_rates[1798];
  //sp 45
  sp_rates[44] += fwd_rates[1798];
  //sp 47
  sp_rates[46] -= fwd_rates[1798];
  //sp 85
  sp_rates[84] += fwd_rates[1798];

  //rxn 1799
  //sp 67
  sp_rates[66] -= fwd_rates[1799];
  //sp 4
  sp_rates[3] += fwd_rates[1799];
  //sp 79
  sp_rates[78] += fwd_rates[1799];
  //sp 80
  sp_rates[79] -= fwd_rates[1799];
  //sp 85
  sp_rates[84] += fwd_rates[1799];

  //rxn 1800
  //sp 67
  sp_rates[66] -= fwd_rates[1800];
  //sp 4
  sp_rates[3] += fwd_rates[1800];
  //sp 48
  sp_rates[47] += fwd_rates[1800];
  //sp 49
  sp_rates[48] -= fwd_rates[1800];
  //sp 85
  sp_rates[84] += fwd_rates[1800];

  //rxn 1801
  //sp 105
  sp_rates[104] += fwd_rates[1801];
  //sp 34
  sp_rates[33] += fwd_rates[1801];
  //sp 35
  sp_rates[34] -= fwd_rates[1801];
  //sp 104
  sp_rates[103] -= fwd_rates[1801];

  //rxn 1802
  //sp 105
  sp_rates[104] += fwd_rates[1802];
  //sp 30
  sp_rates[29] += fwd_rates[1802];
  //sp 31
  sp_rates[30] -= fwd_rates[1802];
  //sp 104
  sp_rates[103] -= fwd_rates[1802];

  //rxn 1803
  //sp 105
  sp_rates[104] += fwd_rates[1803];
  //sp 59
  sp_rates[58] += fwd_rates[1803];
  //sp 61
  sp_rates[60] -= fwd_rates[1803];
  //sp 104
  sp_rates[103] -= fwd_rates[1803];

  //rxn 1804
  //sp 105
  sp_rates[104] += fwd_rates[1804];
  //sp 59
  sp_rates[58] += fwd_rates[1804];
  //sp 60
  sp_rates[59] -= fwd_rates[1804];
  //sp 104
  sp_rates[103] -= fwd_rates[1804];

  //rxn 1805
  //sp 105
  sp_rates[104] += fwd_rates[1805];
  //sp 104
  sp_rates[103] -= fwd_rates[1805];
  //sp 62
  sp_rates[61] += fwd_rates[1805];
  //sp 64
  sp_rates[63] -= fwd_rates[1805];

  //rxn 1806
  //sp 65
  sp_rates[64] -= fwd_rates[1806];
  //sp 105
  sp_rates[104] += fwd_rates[1806];
  //sp 62
  sp_rates[61] += fwd_rates[1806];
  //sp 104
  sp_rates[103] -= fwd_rates[1806];

  //rxn 1807
  //sp 105
  sp_rates[104] += fwd_rates[1807];
  //sp 62
  sp_rates[61] += fwd_rates[1807];
  //sp 63
  sp_rates[62] -= fwd_rates[1807];
  //sp 104
  sp_rates[103] -= fwd_rates[1807];

  //rxn 1808
  //sp 105
  sp_rates[104] += fwd_rates[1808];
  //sp 102
  sp_rates[101] += fwd_rates[1808];
  //sp 103
  sp_rates[102] -= fwd_rates[1808];
  //sp 104
  sp_rates[103] -= fwd_rates[1808];

  //rxn 1809
  //sp 97
  sp_rates[96] += fwd_rates[1809];
  //sp 100
  sp_rates[99] -= fwd_rates[1809];
  //sp 105
  sp_rates[104] += fwd_rates[1809];
  //sp 104
  sp_rates[103] -= fwd_rates[1809];

  //rxn 1810
  //sp 105
  sp_rates[104] += fwd_rates[1810];
  //sp 98
  sp_rates[97] += fwd_rates[1810];
  //sp 99
  sp_rates[98] -= fwd_rates[1810];
  //sp 104
  sp_rates[103] -= fwd_rates[1810];

  //rxn 1811
  //sp 105
  sp_rates[104] += fwd_rates[1811];
  //sp 94
  sp_rates[93] += fwd_rates[1811];
  //sp 95
  sp_rates[94] -= fwd_rates[1811];
  //sp 104
  sp_rates[103] -= fwd_rates[1811];

  //rxn 1812
  //sp 105
  sp_rates[104] += fwd_rates[1812];
  //sp 107
  sp_rates[106] += fwd_rates[1812];
  //sp 109
  sp_rates[108] -= fwd_rates[1812];
  //sp 104
  sp_rates[103] -= fwd_rates[1812];

  //rxn 1813
  //sp 105
  sp_rates[104] += fwd_rates[1813];
  //sp 110
  sp_rates[109] += fwd_rates[1813];
  //sp 111
  sp_rates[110] -= fwd_rates[1813];
  //sp 104
  sp_rates[103] -= fwd_rates[1813];

  //rxn 1814
  //sp 105
  sp_rates[104] += fwd_rates[1814];
  //sp 8
  sp_rates[7] += fwd_rates[1814];
  //sp 6
  sp_rates[5] -= fwd_rates[1814];
  //sp 104
  sp_rates[103] -= fwd_rates[1814];

  //rxn 1815
  //sp 9
  sp_rates[8] += fwd_rates[1815];
  //sp 10
  sp_rates[9] -= fwd_rates[1815];
  //sp 105
  sp_rates[104] += fwd_rates[1815];
  //sp 104
  sp_rates[103] -= fwd_rates[1815];

  //rxn 1816
  //sp 105
  sp_rates[104] += fwd_rates[1816];
  //sp 22
  sp_rates[21] += fwd_rates[1816];
  //sp 23
  sp_rates[22] -= fwd_rates[1816];
  //sp 104
  sp_rates[103] -= fwd_rates[1816];

  //rxn 1817
  //sp 105
  sp_rates[104] += fwd_rates[1817];
  //sp 104
  sp_rates[103] -= fwd_rates[1817];
  //sp 22
  sp_rates[21] += fwd_rates[1817];
  //sp 24
  sp_rates[23] -= fwd_rates[1817];

  //rxn 1818
  //sp 25
  sp_rates[24] += fwd_rates[1818];
  //sp 26
  sp_rates[25] -= fwd_rates[1818];
  //sp 105
  sp_rates[104] += fwd_rates[1818];
  //sp 104
  sp_rates[103] -= fwd_rates[1818];

  //rxn 1819
  //sp 105
  sp_rates[104] += fwd_rates[1819];
  //sp 45
  sp_rates[44] += fwd_rates[1819];
  //sp 46
  sp_rates[45] -= fwd_rates[1819];
  //sp 104
  sp_rates[103] -= fwd_rates[1819];

  //rxn 1820
  //sp 105
  sp_rates[104] += fwd_rates[1820];
  //sp 45
  sp_rates[44] += fwd_rates[1820];
  //sp 47
  sp_rates[46] -= fwd_rates[1820];
  //sp 104
  sp_rates[103] -= fwd_rates[1820];

  //rxn 1821
  //sp 105
  sp_rates[104] += fwd_rates[1821];
  //sp 104
  sp_rates[103] -= fwd_rates[1821];
  //sp 79
  sp_rates[78] += fwd_rates[1821];
  //sp 80
  sp_rates[79] -= fwd_rates[1821];

  //rxn 1822
  //sp 49
  sp_rates[48] -= fwd_rates[1822];
  //sp 105
  sp_rates[104] += fwd_rates[1822];
  //sp 48
  sp_rates[47] += fwd_rates[1822];
  //sp 104
  sp_rates[103] -= fwd_rates[1822];

  //rxn 1823
  //sp 3
  sp_rates[2] += fwd_rates[1823];
  //sp 4
  sp_rates[3] -= fwd_rates[1823];
  //sp 102
  sp_rates[101] -= fwd_rates[1823];
  //sp 103
  sp_rates[102] += fwd_rates[1823];

  //rxn 1824
  //sp 14
  sp_rates[13] += fwd_rates[1824];
  //sp 103
  sp_rates[102] += fwd_rates[1824];
  //sp 102
  sp_rates[101] -= fwd_rates[1824];
  //sp 15
  sp_rates[14] -= fwd_rates[1824];

  //rxn 1825
  //sp 30
  sp_rates[29] += fwd_rates[1825];
  //sp 103
  sp_rates[102] += fwd_rates[1825];
  //sp 102
  sp_rates[101] -= fwd_rates[1825];
  //sp 31
  sp_rates[30] -= fwd_rates[1825];

  //rxn 1826
  //sp 59
  sp_rates[58] += fwd_rates[1826];
  //sp 61
  sp_rates[60] -= fwd_rates[1826];
  //sp 102
  sp_rates[101] -= fwd_rates[1826];
  //sp 103
  sp_rates[102] += fwd_rates[1826];

  //rxn 1827
  //sp 59
  sp_rates[58] += fwd_rates[1827];
  //sp 60
  sp_rates[59] -= fwd_rates[1827];
  //sp 102
  sp_rates[101] -= fwd_rates[1827];
  //sp 103
  sp_rates[102] += fwd_rates[1827];

  //rxn 1828
  //sp 62
  sp_rates[61] += fwd_rates[1828];
  //sp 102
  sp_rates[101] -= fwd_rates[1828];
  //sp 103
  sp_rates[102] += fwd_rates[1828];
  //sp 64
  sp_rates[63] -= fwd_rates[1828];

  //rxn 1829
  //sp 65
  sp_rates[64] -= fwd_rates[1829];
  //sp 62
  sp_rates[61] += fwd_rates[1829];
  //sp 102
  sp_rates[101] -= fwd_rates[1829];
  //sp 103
  sp_rates[102] += fwd_rates[1829];

  //rxn 1830
  //sp 97
  sp_rates[96] += fwd_rates[1830];
  //sp 100
  sp_rates[99] -= fwd_rates[1830];
  //sp 102
  sp_rates[101] -= fwd_rates[1830];
  //sp 103
  sp_rates[102] += fwd_rates[1830];

  //rxn 1831
  //sp 98
  sp_rates[97] += fwd_rates[1831];
  //sp 99
  sp_rates[98] -= fwd_rates[1831];
  //sp 102
  sp_rates[101] -= fwd_rates[1831];
  //sp 103
  sp_rates[102] += fwd_rates[1831];

  //rxn 1832
  //sp 103
  sp_rates[102] += fwd_rates[1832];
  //sp 94
  sp_rates[93] += fwd_rates[1832];
  //sp 102
  sp_rates[101] -= fwd_rates[1832];
  //sp 95
  sp_rates[94] -= fwd_rates[1832];

  //rxn 1833
  //sp 105
  sp_rates[104] -= fwd_rates[1833];
  //sp 102
  sp_rates[101] -= fwd_rates[1833];
  //sp 103
  sp_rates[102] += fwd_rates[1833];
  //sp 104
  sp_rates[103] += fwd_rates[1833];

  //rxn 1834
  //sp 107
  sp_rates[106] += fwd_rates[1834];
  //sp 109
  sp_rates[108] -= fwd_rates[1834];
  //sp 102
  sp_rates[101] -= fwd_rates[1834];
  //sp 103
  sp_rates[102] += fwd_rates[1834];

  //rxn 1835
  //sp 103
  sp_rates[102] += fwd_rates[1835];
  //sp 110
  sp_rates[109] += fwd_rates[1835];
  //sp 102
  sp_rates[101] -= fwd_rates[1835];
  //sp 111
  sp_rates[110] -= fwd_rates[1835];

  //rxn 1836
  //sp 22
  sp_rates[21] += fwd_rates[1836];
  //sp 103
  sp_rates[102] += fwd_rates[1836];
  //sp 102
  sp_rates[101] -= fwd_rates[1836];
  //sp 23
  sp_rates[22] -= fwd_rates[1836];

  //rxn 1837
  //sp 22
  sp_rates[21] += fwd_rates[1837];
  //sp 102
  sp_rates[101] -= fwd_rates[1837];
  //sp 103
  sp_rates[102] += fwd_rates[1837];
  //sp 24
  sp_rates[23] -= fwd_rates[1837];

  //rxn 1838
  //sp 25
  sp_rates[24] += fwd_rates[1838];
  //sp 26
  sp_rates[25] -= fwd_rates[1838];
  //sp 102
  sp_rates[101] -= fwd_rates[1838];
  //sp 103
  sp_rates[102] += fwd_rates[1838];

  //rxn 1839
  //sp 102
  sp_rates[101] -= fwd_rates[1839];
  //sp 45
  sp_rates[44] += fwd_rates[1839];
  //sp 46
  sp_rates[45] -= fwd_rates[1839];
  //sp 103
  sp_rates[102] += fwd_rates[1839];

  //rxn 1840
  //sp 103
  sp_rates[102] += fwd_rates[1840];
  //sp 45
  sp_rates[44] += fwd_rates[1840];
  //sp 102
  sp_rates[101] -= fwd_rates[1840];
  //sp 47
  sp_rates[46] -= fwd_rates[1840];

  //rxn 1841
  //sp 49
  sp_rates[48] -= fwd_rates[1841];
  //sp 102
  sp_rates[101] -= fwd_rates[1841];
  //sp 103
  sp_rates[102] += fwd_rates[1841];
  //sp 48
  sp_rates[47] += fwd_rates[1841];

  //rxn 1842
  //sp 97
  sp_rates[96] -= fwd_rates[1842];
  //sp 3
  sp_rates[2] += fwd_rates[1842];
  //sp 4
  sp_rates[3] -= fwd_rates[1842];
  //sp 99
  sp_rates[98] += fwd_rates[1842];

  //rxn 1843
  //sp 97
  sp_rates[96] -= fwd_rates[1843];
  //sp 99
  sp_rates[98] += fwd_rates[1843];
  //sp 14
  sp_rates[13] += fwd_rates[1843];
  //sp 15
  sp_rates[14] -= fwd_rates[1843];

  //rxn 1844
  //sp 97
  sp_rates[96] -= fwd_rates[1844];
  //sp 34
  sp_rates[33] += fwd_rates[1844];
  //sp 35
  sp_rates[34] -= fwd_rates[1844];
  //sp 99
  sp_rates[98] += fwd_rates[1844];

  //rxn 1845
  //sp 97
  sp_rates[96] -= fwd_rates[1845];
  //sp 99
  sp_rates[98] += fwd_rates[1845];
  //sp 30
  sp_rates[29] += fwd_rates[1845];
  //sp 31
  sp_rates[30] -= fwd_rates[1845];

  //rxn 1846
  //sp 97
  sp_rates[96] -= fwd_rates[1846];
  //sp 59
  sp_rates[58] += fwd_rates[1846];
  //sp 61
  sp_rates[60] -= fwd_rates[1846];
  //sp 99
  sp_rates[98] += fwd_rates[1846];

  //rxn 1847
  //sp 97
  sp_rates[96] -= fwd_rates[1847];
  //sp 59
  sp_rates[58] += fwd_rates[1847];
  //sp 60
  sp_rates[59] -= fwd_rates[1847];
  //sp 99
  sp_rates[98] += fwd_rates[1847];

  //rxn 1848
  //sp 97
  sp_rates[96] -= fwd_rates[1848];
  //sp 99
  sp_rates[98] += fwd_rates[1848];
  //sp 62
  sp_rates[61] += fwd_rates[1848];
  //sp 64
  sp_rates[63] -= fwd_rates[1848];

  //rxn 1849
  //sp 65
  sp_rates[64] -= fwd_rates[1849];
  //sp 99
  sp_rates[98] += fwd_rates[1849];
  //sp 97
  sp_rates[96] -= fwd_rates[1849];
  //sp 62
  sp_rates[61] += fwd_rates[1849];

  //rxn 1850
  //sp 97
  sp_rates[96] -= fwd_rates[1850];
  //sp 99
  sp_rates[98] += fwd_rates[1850];
  //sp 62
  sp_rates[61] += fwd_rates[1850];
  //sp 63
  sp_rates[62] -= fwd_rates[1850];

  //rxn 1851
  //sp 97
  sp_rates[96] -= fwd_rates[1851];
  //sp 99
  sp_rates[98] += fwd_rates[1851];
  //sp 102
  sp_rates[101] += fwd_rates[1851];
  //sp 103
  sp_rates[102] -= fwd_rates[1851];

  //rxn 1852
  //sp 99
  sp_rates[98] += fwd_rates[1852];
  //sp 100
  sp_rates[99] -= fwd_rates[1852];

  //rxn 1853
  //sp 97
  sp_rates[96] -= fwd_rates[1853];
  //sp 98
  sp_rates[97] += fwd_rates[1853];

  //rxn 1854
  //sp 97
  sp_rates[96] -= fwd_rates[1854];
  //sp 99
  sp_rates[98] += fwd_rates[1854];
  //sp 94
  sp_rates[93] += fwd_rates[1854];
  //sp 95
  sp_rates[94] -= fwd_rates[1854];

  //rxn 1855
  //sp 97
  sp_rates[96] -= fwd_rates[1855];
  //sp 99
  sp_rates[98] += fwd_rates[1855];
  //sp 105
  sp_rates[104] -= fwd_rates[1855];
  //sp 104
  sp_rates[103] += fwd_rates[1855];

  //rxn 1856
  //sp 97
  sp_rates[96] -= fwd_rates[1856];
  //sp 99
  sp_rates[98] += fwd_rates[1856];
  //sp 109
  sp_rates[108] -= fwd_rates[1856];
  //sp 107
  sp_rates[106] += fwd_rates[1856];

  //rxn 1857
  //sp 97
  sp_rates[96] -= fwd_rates[1857];
  //sp 99
  sp_rates[98] += fwd_rates[1857];
  //sp 110
  sp_rates[109] += fwd_rates[1857];
  //sp 111
  sp_rates[110] -= fwd_rates[1857];

  //rxn 1858
  //sp 97
  sp_rates[96] -= fwd_rates[1858];
  //sp 99
  sp_rates[98] += fwd_rates[1858];
  //sp 22
  sp_rates[21] += fwd_rates[1858];
  //sp 24
  sp_rates[23] -= fwd_rates[1858];

  //rxn 1859
  //sp 97
  sp_rates[96] -= fwd_rates[1859];
  //sp 26
  sp_rates[25] -= fwd_rates[1859];
  //sp 99
  sp_rates[98] += fwd_rates[1859];
  //sp 25
  sp_rates[24] += fwd_rates[1859];

  //rxn 1860
  //sp 97
  sp_rates[96] -= fwd_rates[1860];
  //sp 99
  sp_rates[98] += fwd_rates[1860];
  //sp 45
  sp_rates[44] += fwd_rates[1860];
  //sp 46
  sp_rates[45] -= fwd_rates[1860];

  //rxn 1861
  //sp 97
  sp_rates[96] -= fwd_rates[1861];
  //sp 99
  sp_rates[98] += fwd_rates[1861];
  //sp 45
  sp_rates[44] += fwd_rates[1861];
  //sp 47
  sp_rates[46] -= fwd_rates[1861];

  //rxn 1862
  //sp 97
  sp_rates[96] -= fwd_rates[1862];
  //sp 99
  sp_rates[98] += fwd_rates[1862];
  //sp 79
  sp_rates[78] += fwd_rates[1862];
  //sp 80
  sp_rates[79] -= fwd_rates[1862];

  //rxn 1863
  //sp 49
  sp_rates[48] -= fwd_rates[1863];
  //sp 99
  sp_rates[98] += fwd_rates[1863];
  //sp 97
  sp_rates[96] -= fwd_rates[1863];
  //sp 48
  sp_rates[47] += fwd_rates[1863];

  //rxn 1864
  //sp 97
  sp_rates[96] -= fwd_rates[1864];
  //sp 100
  sp_rates[99] += fwd_rates[1864];
  //sp 3
  sp_rates[2] += fwd_rates[1864];
  //sp 4
  sp_rates[3] -= fwd_rates[1864];

  //rxn 1865
  //sp 97
  sp_rates[96] -= fwd_rates[1865];
  //sp 100
  sp_rates[99] += fwd_rates[1865];
  //sp 14
  sp_rates[13] += fwd_rates[1865];
  //sp 15
  sp_rates[14] -= fwd_rates[1865];

  //rxn 1866
  //sp 97
  sp_rates[96] -= fwd_rates[1866];
  //sp 34
  sp_rates[33] += fwd_rates[1866];
  //sp 35
  sp_rates[34] -= fwd_rates[1866];
  //sp 100
  sp_rates[99] += fwd_rates[1866];

  //rxn 1867
  //sp 97
  sp_rates[96] -= fwd_rates[1867];
  //sp 100
  sp_rates[99] += fwd_rates[1867];
  //sp 30
  sp_rates[29] += fwd_rates[1867];
  //sp 31
  sp_rates[30] -= fwd_rates[1867];

  //rxn 1868
  //sp 97
  sp_rates[96] -= fwd_rates[1868];
  //sp 59
  sp_rates[58] += fwd_rates[1868];
  //sp 100
  sp_rates[99] += fwd_rates[1868];
  //sp 61
  sp_rates[60] -= fwd_rates[1868];

  //rxn 1869
  //sp 97
  sp_rates[96] -= fwd_rates[1869];
  //sp 100
  sp_rates[99] += fwd_rates[1869];
  //sp 59
  sp_rates[58] += fwd_rates[1869];
  //sp 60
  sp_rates[59] -= fwd_rates[1869];

  //rxn 1870
  //sp 97
  sp_rates[96] -= fwd_rates[1870];
  //sp 100
  sp_rates[99] += fwd_rates[1870];
  //sp 62
  sp_rates[61] += fwd_rates[1870];
  //sp 64
  sp_rates[63] -= fwd_rates[1870];

  //rxn 1871
  //sp 65
  sp_rates[64] -= fwd_rates[1871];
  //sp 100
  sp_rates[99] += fwd_rates[1871];
  //sp 97
  sp_rates[96] -= fwd_rates[1871];
  //sp 62
  sp_rates[61] += fwd_rates[1871];

  //rxn 1872
  //sp 97
  sp_rates[96] -= fwd_rates[1872];
  //sp 100
  sp_rates[99] += fwd_rates[1872];
  //sp 62
  sp_rates[61] += fwd_rates[1872];
  //sp 63
  sp_rates[62] -= fwd_rates[1872];

  //rxn 1873
  //sp 97
  sp_rates[96] -= fwd_rates[1873];
  //sp 100
  sp_rates[99] += fwd_rates[1873];
  //sp 102
  sp_rates[101] += fwd_rates[1873];
  //sp 103
  sp_rates[102] -= fwd_rates[1873];

  //rxn 1874
  //sp 97
  sp_rates[96] -= fwd_rates[1874];
  //sp 98
  sp_rates[97] += fwd_rates[1874];
  //sp 99
  sp_rates[98] -= fwd_rates[1874];
  //sp 100
  sp_rates[99] += fwd_rates[1874];

  //rxn 1875
  //sp 97
  sp_rates[96] -= fwd_rates[1875];
  //sp 100
  sp_rates[99] += fwd_rates[1875];
  //sp 94
  sp_rates[93] += fwd_rates[1875];
  //sp 95
  sp_rates[94] -= fwd_rates[1875];

  //rxn 1876
  //sp 97
  sp_rates[96] -= fwd_rates[1876];
  //sp 100
  sp_rates[99] += fwd_rates[1876];
  //sp 105
  sp_rates[104] -= fwd_rates[1876];
  //sp 104
  sp_rates[103] += fwd_rates[1876];

  //rxn 1877
  //sp 97
  sp_rates[96] -= fwd_rates[1877];
  //sp 107
  sp_rates[106] += fwd_rates[1877];
  //sp 100
  sp_rates[99] += fwd_rates[1877];
  //sp 109
  sp_rates[108] -= fwd_rates[1877];

  //rxn 1878
  //sp 97
  sp_rates[96] -= fwd_rates[1878];
  //sp 100
  sp_rates[99] += fwd_rates[1878];
  //sp 110
  sp_rates[109] += fwd_rates[1878];
  //sp 111
  sp_rates[110] -= fwd_rates[1878];

  //rxn 1879
  //sp 97
  sp_rates[96] -= fwd_rates[1879];
  //sp 100
  sp_rates[99] += fwd_rates[1879];
  //sp 22
  sp_rates[21] += fwd_rates[1879];
  //sp 24
  sp_rates[23] -= fwd_rates[1879];

  //rxn 1880
  //sp 97
  sp_rates[96] -= fwd_rates[1880];
  //sp 26
  sp_rates[25] -= fwd_rates[1880];
  //sp 100
  sp_rates[99] += fwd_rates[1880];
  //sp 25
  sp_rates[24] += fwd_rates[1880];

  //rxn 1881
  //sp 97
  sp_rates[96] -= fwd_rates[1881];
  //sp 100
  sp_rates[99] += fwd_rates[1881];
  //sp 45
  sp_rates[44] += fwd_rates[1881];
  //sp 46
  sp_rates[45] -= fwd_rates[1881];

  //rxn 1882
  //sp 97
  sp_rates[96] -= fwd_rates[1882];
  //sp 100
  sp_rates[99] += fwd_rates[1882];
  //sp 45
  sp_rates[44] += fwd_rates[1882];
  //sp 47
  sp_rates[46] -= fwd_rates[1882];

  //rxn 1883
  //sp 97
  sp_rates[96] -= fwd_rates[1883];
  //sp 100
  sp_rates[99] += fwd_rates[1883];
  //sp 79
  sp_rates[78] += fwd_rates[1883];
  //sp 80
  sp_rates[79] -= fwd_rates[1883];

  //rxn 1884
  //sp 49
  sp_rates[48] -= fwd_rates[1884];
  //sp 100
  sp_rates[99] += fwd_rates[1884];
  //sp 97
  sp_rates[96] -= fwd_rates[1884];
  //sp 48
  sp_rates[47] += fwd_rates[1884];

  //rxn 1885
  //sp 98
  sp_rates[97] -= fwd_rates[1885];
  //sp 3
  sp_rates[2] += fwd_rates[1885];
  //sp 4
  sp_rates[3] -= fwd_rates[1885];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1885];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1885];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1885];

  //rxn 1886
  //sp 98
  sp_rates[97] -= fwd_rates[1886];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1886];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1886];
  //sp 14
  sp_rates[13] += fwd_rates[1886];
  //sp 15
  sp_rates[14] -= 0.97 * fwd_rates[1886];

  //rxn 1887
  //sp 98
  sp_rates[97] -= fwd_rates[1887];
  //sp 35
  sp_rates[34] -= fwd_rates[1887];
  //sp 34
  sp_rates[33] += fwd_rates[1887];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1887];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1887];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1887];

  //rxn 1888
  //sp 98
  sp_rates[97] -= fwd_rates[1888];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1888];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1888];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1888];
  //sp 30
  sp_rates[29] += fwd_rates[1888];
  //sp 31
  sp_rates[30] -= fwd_rates[1888];

  //rxn 1889
  //sp 98
  sp_rates[97] -= fwd_rates[1889];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1889];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1889];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1889];
  //sp 59
  sp_rates[58] += fwd_rates[1889];
  //sp 61
  sp_rates[60] -= fwd_rates[1889];

  //rxn 1890
  //sp 98
  sp_rates[97] -= fwd_rates[1890];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1890];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1890];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1890];
  //sp 59
  sp_rates[58] += fwd_rates[1890];
  //sp 60
  sp_rates[59] -= fwd_rates[1890];

  //rxn 1891
  //sp 98
  sp_rates[97] -= fwd_rates[1891];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1891];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1891];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1891];
  //sp 62
  sp_rates[61] += fwd_rates[1891];
  //sp 64
  sp_rates[63] -= fwd_rates[1891];

  //rxn 1892
  //sp 65
  sp_rates[64] -= fwd_rates[1892];
  //sp 98
  sp_rates[97] -= fwd_rates[1892];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1892];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1892];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1892];
  //sp 62
  sp_rates[61] += fwd_rates[1892];

  //rxn 1893
  //sp 98
  sp_rates[97] -= fwd_rates[1893];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1893];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1893];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1893];
  //sp 62
  sp_rates[61] += fwd_rates[1893];
  //sp 63
  sp_rates[62] -= fwd_rates[1893];

  //rxn 1894
  //sp 98
  sp_rates[97] -= fwd_rates[1894];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1894];
  //sp 102
  sp_rates[101] += fwd_rates[1894];
  //sp 103
  sp_rates[102] -= fwd_rates[1894];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1894];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1894];

  //rxn 1895
  //sp 97
  sp_rates[96] += fwd_rates[1895];
  //sp 98
  sp_rates[97] -= fwd_rates[1895];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1895];
  //sp 100
  sp_rates[99] -= fwd_rates[1895];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1895];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1895];

  //rxn 1896
  //sp 99
  sp_rates[98] -= 0.030000000000000027 * fwd_rates[1896];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1896];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1896];

  //rxn 1897
  //sp 98
  sp_rates[97] -= fwd_rates[1897];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1897];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1897];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1897];
  //sp 94
  sp_rates[93] += fwd_rates[1897];
  //sp 95
  sp_rates[94] -= fwd_rates[1897];

  //rxn 1898
  //sp 98
  sp_rates[97] -= fwd_rates[1898];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1898];
  //sp 104
  sp_rates[103] += fwd_rates[1898];
  //sp 105
  sp_rates[104] -= fwd_rates[1898];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1898];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1898];

  //rxn 1899
  //sp 98
  sp_rates[97] -= fwd_rates[1899];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1899];
  //sp 107
  sp_rates[106] += fwd_rates[1899];
  //sp 109
  sp_rates[108] -= fwd_rates[1899];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1899];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1899];

  //rxn 1900
  //sp 98
  sp_rates[97] -= fwd_rates[1900];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1900];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1900];
  //sp 111
  sp_rates[110] -= fwd_rates[1900];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1900];
  //sp 110
  sp_rates[109] += fwd_rates[1900];

  //rxn 1901
  //sp 98
  sp_rates[97] -= fwd_rates[1901];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1901];
  //sp 6
  sp_rates[5] -= fwd_rates[1901];
  //sp 8
  sp_rates[7] += fwd_rates[1901];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1901];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1901];

  //rxn 1902
  //sp 98
  sp_rates[97] -= fwd_rates[1902];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1902];
  //sp 5
  sp_rates[4] -= fwd_rates[1902];
  //sp 10
  sp_rates[9] += fwd_rates[1902];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1902];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1902];

  //rxn 1903
  //sp 98
  sp_rates[97] -= fwd_rates[1903];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1903];
  //sp 7
  sp_rates[6] += fwd_rates[1903];
  //sp 8
  sp_rates[7] -= fwd_rates[1903];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1903];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1903];

  //rxn 1904
  //sp 98
  sp_rates[97] -= fwd_rates[1904];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1904];
  //sp 9
  sp_rates[8] += fwd_rates[1904];
  //sp 10
  sp_rates[9] -= fwd_rates[1904];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1904];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1904];

  //rxn 1905
  //sp 98
  sp_rates[97] -= fwd_rates[1905];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1905];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1905];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1905];
  //sp 22
  sp_rates[21] += fwd_rates[1905];
  //sp 23
  sp_rates[22] -= fwd_rates[1905];

  //rxn 1906
  //sp 98
  sp_rates[97] -= fwd_rates[1906];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1906];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1906];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1906];
  //sp 22
  sp_rates[21] += fwd_rates[1906];
  //sp 24
  sp_rates[23] -= fwd_rates[1906];

  //rxn 1907
  //sp 98
  sp_rates[97] -= fwd_rates[1907];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1907];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1907];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1907];
  //sp 25
  sp_rates[24] += fwd_rates[1907];
  //sp 26
  sp_rates[25] -= fwd_rates[1907];

  //rxn 1908
  //sp 98
  sp_rates[97] -= fwd_rates[1908];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1908];
  //sp 45
  sp_rates[44] += fwd_rates[1908];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1908];
  //sp 46
  sp_rates[45] -= fwd_rates[1908];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1908];

  //rxn 1909
  //sp 98
  sp_rates[97] -= fwd_rates[1909];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1909];
  //sp 45
  sp_rates[44] += fwd_rates[1909];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1909];
  //sp 47
  sp_rates[46] -= fwd_rates[1909];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1909];

  //rxn 1910
  //sp 98
  sp_rates[97] -= fwd_rates[1910];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1910];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1910];
  //sp 79
  sp_rates[78] += fwd_rates[1910];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1910];
  //sp 80
  sp_rates[79] -= fwd_rates[1910];

  //rxn 1911
  //sp 98
  sp_rates[97] -= fwd_rates[1911];
  //sp 99
  sp_rates[98] += 0.97 * fwd_rates[1911];
  //sp 78
  sp_rates[77] += 0.03 * fwd_rates[1911];
  //sp 15
  sp_rates[14] += 0.03 * fwd_rates[1911];
  //sp 48
  sp_rates[47] += fwd_rates[1911];
  //sp 49
  sp_rates[48] -= fwd_rates[1911];

  //rxn 1912
  //sp 3
  sp_rates[2] += fwd_rates[1912];
  //sp 4
  sp_rates[3] -= fwd_rates[1912];
  //sp 94
  sp_rates[93] -= fwd_rates[1912];
  //sp 95
  sp_rates[94] += fwd_rates[1912];

  //rxn 1913
  //sp 95
  sp_rates[94] += fwd_rates[1913];
  //sp 14
  sp_rates[13] += fwd_rates[1913];
  //sp 94
  sp_rates[93] -= fwd_rates[1913];
  //sp 15
  sp_rates[14] -= fwd_rates[1913];

  //rxn 1914
  //sp 34
  sp_rates[33] += fwd_rates[1914];
  //sp 35
  sp_rates[34] -= fwd_rates[1914];
  //sp 94
  sp_rates[93] -= fwd_rates[1914];
  //sp 95
  sp_rates[94] += fwd_rates[1914];

  //rxn 1915
  //sp 95
  sp_rates[94] += fwd_rates[1915];
  //sp 30
  sp_rates[29] += fwd_rates[1915];
  //sp 94
  sp_rates[93] -= fwd_rates[1915];
  //sp 31
  sp_rates[30] -= fwd_rates[1915];

  //rxn 1916
  //sp 59
  sp_rates[58] += fwd_rates[1916];
  //sp 61
  sp_rates[60] -= fwd_rates[1916];
  //sp 94
  sp_rates[93] -= fwd_rates[1916];
  //sp 95
  sp_rates[94] += fwd_rates[1916];

  //rxn 1917
  //sp 59
  sp_rates[58] += fwd_rates[1917];
  //sp 60
  sp_rates[59] -= fwd_rates[1917];
  //sp 94
  sp_rates[93] -= fwd_rates[1917];
  //sp 95
  sp_rates[94] += fwd_rates[1917];

  //rxn 1918
  //sp 95
  sp_rates[94] += fwd_rates[1918];
  //sp 102
  sp_rates[101] += fwd_rates[1918];
  //sp 94
  sp_rates[93] -= fwd_rates[1918];
  //sp 103
  sp_rates[102] -= fwd_rates[1918];

  //rxn 1919
  //sp 97
  sp_rates[96] += fwd_rates[1919];
  //sp 100
  sp_rates[99] -= fwd_rates[1919];
  //sp 94
  sp_rates[93] -= fwd_rates[1919];
  //sp 95
  sp_rates[94] += fwd_rates[1919];

  //rxn 1920
  //sp 98
  sp_rates[97] += fwd_rates[1920];
  //sp 99
  sp_rates[98] -= fwd_rates[1920];
  //sp 94
  sp_rates[93] -= fwd_rates[1920];
  //sp 95
  sp_rates[94] += fwd_rates[1920];

  //rxn 1921
  //sp 105
  sp_rates[104] -= fwd_rates[1921];
  //sp 94
  sp_rates[93] -= fwd_rates[1921];
  //sp 95
  sp_rates[94] += fwd_rates[1921];
  //sp 104
  sp_rates[103] += fwd_rates[1921];

  //rxn 1922
  //sp 107
  sp_rates[106] += fwd_rates[1922];
  //sp 109
  sp_rates[108] -= fwd_rates[1922];
  //sp 94
  sp_rates[93] -= fwd_rates[1922];
  //sp 95
  sp_rates[94] += fwd_rates[1922];

  //rxn 1923
  //sp 95
  sp_rates[94] += fwd_rates[1923];
  //sp 110
  sp_rates[109] += fwd_rates[1923];
  //sp 94
  sp_rates[93] -= fwd_rates[1923];
  //sp 111
  sp_rates[110] -= fwd_rates[1923];

  //rxn 1924
  //sp 22
  sp_rates[21] += fwd_rates[1924];
  //sp 94
  sp_rates[93] -= fwd_rates[1924];
  //sp 95
  sp_rates[94] += fwd_rates[1924];
  //sp 24
  sp_rates[23] -= fwd_rates[1924];

  //rxn 1925
  //sp 25
  sp_rates[24] += fwd_rates[1925];
  //sp 26
  sp_rates[25] -= fwd_rates[1925];
  //sp 94
  sp_rates[93] -= fwd_rates[1925];
  //sp 95
  sp_rates[94] += fwd_rates[1925];

  //rxn 1926
  //sp 95
  sp_rates[94] += fwd_rates[1926];
  //sp 94
  sp_rates[93] -= fwd_rates[1926];
  //sp 46
  sp_rates[45] -= fwd_rates[1926];
  //sp 45
  sp_rates[44] += fwd_rates[1926];

  //rxn 1927
  //sp 95
  sp_rates[94] += fwd_rates[1927];
  //sp 45
  sp_rates[44] += fwd_rates[1927];
  //sp 94
  sp_rates[93] -= fwd_rates[1927];
  //sp 47
  sp_rates[46] -= fwd_rates[1927];

  //rxn 1928
  //sp 95
  sp_rates[94] += fwd_rates[1928];
  //sp 94
  sp_rates[93] -= fwd_rates[1928];
  //sp 79
  sp_rates[78] += fwd_rates[1928];
  //sp 80
  sp_rates[79] -= fwd_rates[1928];

  //rxn 1929
  //sp 49
  sp_rates[48] -= fwd_rates[1929];
  //sp 94
  sp_rates[93] -= fwd_rates[1929];
  //sp 95
  sp_rates[94] += fwd_rates[1929];
  //sp 48
  sp_rates[47] += fwd_rates[1929];

  //rxn 1930
  //sp 30
  sp_rates[29] += fwd_rates[1930];
  //sp 111
  sp_rates[110] += fwd_rates[1930];
  //sp 110
  sp_rates[109] -= fwd_rates[1930];
  //sp 31
  sp_rates[30] -= fwd_rates[1930];

  //rxn 1931
  //sp 59
  sp_rates[58] += fwd_rates[1931];
  //sp 61
  sp_rates[60] -= fwd_rates[1931];
  //sp 110
  sp_rates[109] -= fwd_rates[1931];
  //sp 111
  sp_rates[110] += fwd_rates[1931];

  //rxn 1932
  //sp 59
  sp_rates[58] += fwd_rates[1932];
  //sp 60
  sp_rates[59] -= fwd_rates[1932];
  //sp 110
  sp_rates[109] -= fwd_rates[1932];
  //sp 111
  sp_rates[110] += fwd_rates[1932];

  //rxn 1933
  //sp 62
  sp_rates[61] += fwd_rates[1933];
  //sp 110
  sp_rates[109] -= fwd_rates[1933];
  //sp 111
  sp_rates[110] += fwd_rates[1933];
  //sp 64
  sp_rates[63] -= fwd_rates[1933];

  //rxn 1934
  //sp 65
  sp_rates[64] -= fwd_rates[1934];
  //sp 62
  sp_rates[61] += fwd_rates[1934];
  //sp 110
  sp_rates[109] -= fwd_rates[1934];
  //sp 111
  sp_rates[110] += fwd_rates[1934];

  //rxn 1935
  //sp 111
  sp_rates[110] += fwd_rates[1935];
  //sp 62
  sp_rates[61] += fwd_rates[1935];
  //sp 110
  sp_rates[109] -= fwd_rates[1935];
  //sp 63
  sp_rates[62] -= fwd_rates[1935];

  //rxn 1936
  //sp 102
  sp_rates[101] += fwd_rates[1936];
  //sp 111
  sp_rates[110] += fwd_rates[1936];
  //sp 110
  sp_rates[109] -= fwd_rates[1936];
  //sp 103
  sp_rates[102] -= fwd_rates[1936];

  //rxn 1937
  //sp 97
  sp_rates[96] += fwd_rates[1937];
  //sp 100
  sp_rates[99] -= fwd_rates[1937];
  //sp 110
  sp_rates[109] -= fwd_rates[1937];
  //sp 111
  sp_rates[110] += fwd_rates[1937];

  //rxn 1938
  //sp 98
  sp_rates[97] += fwd_rates[1938];
  //sp 99
  sp_rates[98] -= fwd_rates[1938];
  //sp 110
  sp_rates[109] -= fwd_rates[1938];
  //sp 111
  sp_rates[110] += fwd_rates[1938];

  //rxn 1939
  //sp 111
  sp_rates[110] += fwd_rates[1939];
  //sp 94
  sp_rates[93] += fwd_rates[1939];
  //sp 110
  sp_rates[109] -= fwd_rates[1939];
  //sp 95
  sp_rates[94] -= fwd_rates[1939];

  //rxn 1940
  //sp 105
  sp_rates[104] -= fwd_rates[1940];
  //sp 110
  sp_rates[109] -= fwd_rates[1940];
  //sp 111
  sp_rates[110] += fwd_rates[1940];
  //sp 104
  sp_rates[103] += fwd_rates[1940];

  //rxn 1941
  //sp 107
  sp_rates[106] += fwd_rates[1941];
  //sp 109
  sp_rates[108] -= fwd_rates[1941];
  //sp 110
  sp_rates[109] -= fwd_rates[1941];
  //sp 111
  sp_rates[110] += fwd_rates[1941];

  //rxn 1942
  //sp 110
  sp_rates[109] -= fwd_rates[1942];
  //sp 6
  sp_rates[5] -= fwd_rates[1942];
  //sp 111
  sp_rates[110] += fwd_rates[1942];
  //sp 8
  sp_rates[7] += fwd_rates[1942];

  //rxn 1943
  //sp 10
  sp_rates[9] += fwd_rates[1943];
  //sp 5
  sp_rates[4] -= fwd_rates[1943];
  //sp 110
  sp_rates[109] -= fwd_rates[1943];
  //sp 111
  sp_rates[110] += fwd_rates[1943];

  //rxn 1944
  //sp 9
  sp_rates[8] += fwd_rates[1944];
  //sp 10
  sp_rates[9] -= fwd_rates[1944];
  //sp 110
  sp_rates[109] -= fwd_rates[1944];
  //sp 111
  sp_rates[110] += fwd_rates[1944];

  //rxn 1945
  //sp 22
  sp_rates[21] += fwd_rates[1945];
  //sp 111
  sp_rates[110] += fwd_rates[1945];
  //sp 110
  sp_rates[109] -= fwd_rates[1945];
  //sp 23
  sp_rates[22] -= fwd_rates[1945];

  //rxn 1946
  //sp 22
  sp_rates[21] += fwd_rates[1946];
  //sp 110
  sp_rates[109] -= fwd_rates[1946];
  //sp 111
  sp_rates[110] += fwd_rates[1946];
  //sp 24
  sp_rates[23] -= fwd_rates[1946];

  //rxn 1947
  //sp 25
  sp_rates[24] += fwd_rates[1947];
  //sp 26
  sp_rates[25] -= fwd_rates[1947];
  //sp 110
  sp_rates[109] -= fwd_rates[1947];
  //sp 111
  sp_rates[110] += fwd_rates[1947];

  //rxn 1948
  //sp 110
  sp_rates[109] -= fwd_rates[1948];
  //sp 45
  sp_rates[44] += fwd_rates[1948];
  //sp 46
  sp_rates[45] -= fwd_rates[1948];
  //sp 111
  sp_rates[110] += fwd_rates[1948];

  //rxn 1949
  //sp 111
  sp_rates[110] += fwd_rates[1949];
  //sp 45
  sp_rates[44] += fwd_rates[1949];
  //sp 110
  sp_rates[109] -= fwd_rates[1949];
  //sp 47
  sp_rates[46] -= fwd_rates[1949];

  //rxn 1950
  //sp 111
  sp_rates[110] += fwd_rates[1950];
  //sp 110
  sp_rates[109] -= fwd_rates[1950];
  //sp 79
  sp_rates[78] += fwd_rates[1950];
  //sp 80
  sp_rates[79] -= fwd_rates[1950];

  //rxn 1951
  //sp 49
  sp_rates[48] -= fwd_rates[1951];
  //sp 110
  sp_rates[109] -= fwd_rates[1951];
  //sp 111
  sp_rates[110] += fwd_rates[1951];
  //sp 48
  sp_rates[47] += fwd_rates[1951];

  //rxn 1952
  //sp 107
  sp_rates[106] -= fwd_rates[1952];
  //sp 109
  sp_rates[108] += fwd_rates[1952];
  //sp 14
  sp_rates[13] += fwd_rates[1952];
  //sp 15
  sp_rates[14] -= fwd_rates[1952];

  //rxn 1953
  //sp 34
  sp_rates[33] += fwd_rates[1953];
  //sp 35
  sp_rates[34] -= fwd_rates[1953];
  //sp 109
  sp_rates[108] += fwd_rates[1953];
  //sp 107
  sp_rates[106] -= fwd_rates[1953];

  //rxn 1954
  //sp 107
  sp_rates[106] -= fwd_rates[1954];
  //sp 109
  sp_rates[108] += fwd_rates[1954];
  //sp 30
  sp_rates[29] += fwd_rates[1954];
  //sp 31
  sp_rates[30] -= fwd_rates[1954];

  //rxn 1955
  //sp 109
  sp_rates[108] += fwd_rates[1955];
  //sp 107
  sp_rates[106] -= fwd_rates[1955];
  //sp 61
  sp_rates[60] -= fwd_rates[1955];
  //sp 59
  sp_rates[58] += fwd_rates[1955];

  //rxn 1956
  //sp 109
  sp_rates[108] += fwd_rates[1956];
  //sp 107
  sp_rates[106] -= fwd_rates[1956];
  //sp 60
  sp_rates[59] -= fwd_rates[1956];
  //sp 59
  sp_rates[58] += fwd_rates[1956];

  //rxn 1957
  //sp 107
  sp_rates[106] -= fwd_rates[1957];
  //sp 109
  sp_rates[108] += fwd_rates[1957];
  //sp 62
  sp_rates[61] += fwd_rates[1957];
  //sp 64
  sp_rates[63] -= fwd_rates[1957];

  //rxn 1958
  //sp 65
  sp_rates[64] -= fwd_rates[1958];
  //sp 107
  sp_rates[106] -= fwd_rates[1958];
  //sp 109
  sp_rates[108] += fwd_rates[1958];
  //sp 62
  sp_rates[61] += fwd_rates[1958];

  //rxn 1959
  //sp 107
  sp_rates[106] -= fwd_rates[1959];
  //sp 109
  sp_rates[108] += fwd_rates[1959];
  //sp 62
  sp_rates[61] += fwd_rates[1959];
  //sp 63
  sp_rates[62] -= fwd_rates[1959];

  //rxn 1960
  //sp 107
  sp_rates[106] -= fwd_rates[1960];
  //sp 109
  sp_rates[108] += fwd_rates[1960];
  //sp 102
  sp_rates[101] += fwd_rates[1960];
  //sp 103
  sp_rates[102] -= fwd_rates[1960];

  //rxn 1961
  //sp 97
  sp_rates[96] += fwd_rates[1961];
  //sp 107
  sp_rates[106] -= fwd_rates[1961];
  //sp 100
  sp_rates[99] -= fwd_rates[1961];
  //sp 109
  sp_rates[108] += fwd_rates[1961];

  //rxn 1962
  //sp 98
  sp_rates[97] += fwd_rates[1962];
  //sp 99
  sp_rates[98] -= fwd_rates[1962];
  //sp 109
  sp_rates[108] += fwd_rates[1962];
  //sp 107
  sp_rates[106] -= fwd_rates[1962];

  //rxn 1963
  //sp 107
  sp_rates[106] -= fwd_rates[1963];
  //sp 109
  sp_rates[108] += fwd_rates[1963];
  //sp 94
  sp_rates[93] += fwd_rates[1963];
  //sp 95
  sp_rates[94] -= fwd_rates[1963];

  //rxn 1964
  //sp 105
  sp_rates[104] -= fwd_rates[1964];
  //sp 107
  sp_rates[106] -= fwd_rates[1964];
  //sp 109
  sp_rates[108] += fwd_rates[1964];
  //sp 104
  sp_rates[103] += fwd_rates[1964];

  //rxn 1965
  //sp 107
  sp_rates[106] -= fwd_rates[1965];
  //sp 109
  sp_rates[108] += fwd_rates[1965];
  //sp 110
  sp_rates[109] += fwd_rates[1965];
  //sp 111
  sp_rates[110] -= fwd_rates[1965];

  //rxn 1966
  //sp 107
  sp_rates[106] -= fwd_rates[1966];
  //sp 109
  sp_rates[108] += fwd_rates[1966];
  //sp 22
  sp_rates[21] += fwd_rates[1966];
  //sp 23
  sp_rates[22] -= fwd_rates[1966];

  //rxn 1967
  //sp 107
  sp_rates[106] -= fwd_rates[1967];
  //sp 109
  sp_rates[108] += fwd_rates[1967];
  //sp 22
  sp_rates[21] += fwd_rates[1967];
  //sp 24
  sp_rates[23] -= fwd_rates[1967];

  //rxn 1968
  //sp 109
  sp_rates[108] += fwd_rates[1968];
  //sp 107
  sp_rates[106] -= fwd_rates[1968];
  //sp 45
  sp_rates[44] += fwd_rates[1968];
  //sp 46
  sp_rates[45] -= fwd_rates[1968];

  //rxn 1969
  //sp 109
  sp_rates[108] += fwd_rates[1969];
  //sp 107
  sp_rates[106] -= fwd_rates[1969];
  //sp 45
  sp_rates[44] += fwd_rates[1969];
  //sp 47
  sp_rates[46] -= fwd_rates[1969];

  //rxn 1970
  //sp 107
  sp_rates[106] -= fwd_rates[1970];
  //sp 109
  sp_rates[108] += fwd_rates[1970];
  //sp 79
  sp_rates[78] += fwd_rates[1970];
  //sp 80
  sp_rates[79] -= fwd_rates[1970];

  //rxn 1971
  //sp 49
  sp_rates[48] -= fwd_rates[1971];
  //sp 107
  sp_rates[106] -= fwd_rates[1971];
  //sp 109
  sp_rates[108] += fwd_rates[1971];
  //sp 48
  sp_rates[47] += fwd_rates[1971];

  //rxn 1972
  //sp 3
  sp_rates[2] += fwd_rates[1972];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1972];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1972];
  //sp 112
  sp_rates[111] -= fwd_rates[1972];

  //rxn 1973
  //sp 34
  sp_rates[33] += fwd_rates[1973];
  //sp 35
  sp_rates[34] -= fwd_rates[1973];
  //sp 4
  sp_rates[3] += fwd_rates[1973];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1973];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1973];
  //sp 112
  sp_rates[111] -= fwd_rates[1973];

  //rxn 1974
  //sp 4
  sp_rates[3] += fwd_rates[1974];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1974];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1974];
  //sp 112
  sp_rates[111] -= fwd_rates[1974];
  //sp 30
  sp_rates[29] += fwd_rates[1974];
  //sp 31
  sp_rates[30] -= fwd_rates[1974];

  //rxn 1975
  //sp 4
  sp_rates[3] += fwd_rates[1975];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1975];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1975];
  //sp 112
  sp_rates[111] -= fwd_rates[1975];
  //sp 59
  sp_rates[58] += fwd_rates[1975];
  //sp 61
  sp_rates[60] -= fwd_rates[1975];

  //rxn 1976
  //sp 4
  sp_rates[3] += fwd_rates[1976];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1976];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1976];
  //sp 112
  sp_rates[111] -= fwd_rates[1976];
  //sp 59
  sp_rates[58] += fwd_rates[1976];
  //sp 60
  sp_rates[59] -= fwd_rates[1976];

  //rxn 1977
  //sp 4
  sp_rates[3] += fwd_rates[1977];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1977];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1977];
  //sp 112
  sp_rates[111] -= fwd_rates[1977];
  //sp 62
  sp_rates[61] += fwd_rates[1977];
  //sp 64
  sp_rates[63] -= fwd_rates[1977];

  //rxn 1978
  //sp 65
  sp_rates[64] -= fwd_rates[1978];
  //sp 4
  sp_rates[3] += fwd_rates[1978];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1978];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1978];
  //sp 112
  sp_rates[111] -= fwd_rates[1978];
  //sp 62
  sp_rates[61] += fwd_rates[1978];

  //rxn 1979
  //sp 4
  sp_rates[3] += fwd_rates[1979];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1979];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1979];
  //sp 112
  sp_rates[111] -= fwd_rates[1979];
  //sp 62
  sp_rates[61] += fwd_rates[1979];
  //sp 63
  sp_rates[62] -= fwd_rates[1979];

  //rxn 1980
  //sp 4
  sp_rates[3] += fwd_rates[1980];
  //sp 102
  sp_rates[101] += fwd_rates[1980];
  //sp 103
  sp_rates[102] -= fwd_rates[1980];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1980];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1980];
  //sp 112
  sp_rates[111] -= fwd_rates[1980];

  //rxn 1981
  //sp 97
  sp_rates[96] += fwd_rates[1981];
  //sp 4
  sp_rates[3] += fwd_rates[1981];
  //sp 100
  sp_rates[99] -= fwd_rates[1981];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1981];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1981];
  //sp 112
  sp_rates[111] -= fwd_rates[1981];

  //rxn 1982
  //sp 98
  sp_rates[97] += fwd_rates[1982];
  //sp 99
  sp_rates[98] -= fwd_rates[1982];
  //sp 4
  sp_rates[3] += fwd_rates[1982];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1982];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1982];
  //sp 112
  sp_rates[111] -= fwd_rates[1982];

  //rxn 1983
  //sp 4
  sp_rates[3] += fwd_rates[1983];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1983];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1983];
  //sp 112
  sp_rates[111] -= fwd_rates[1983];
  //sp 94
  sp_rates[93] += fwd_rates[1983];
  //sp 95
  sp_rates[94] -= fwd_rates[1983];

  //rxn 1984
  //sp 4
  sp_rates[3] += fwd_rates[1984];
  //sp 104
  sp_rates[103] += fwd_rates[1984];
  //sp 105
  sp_rates[104] -= fwd_rates[1984];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1984];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1984];
  //sp 112
  sp_rates[111] -= fwd_rates[1984];

  //rxn 1985
  //sp 4
  sp_rates[3] += fwd_rates[1985];
  //sp 107
  sp_rates[106] += 1.85 * fwd_rates[1985];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1985];
  //sp 109
  sp_rates[108] -= fwd_rates[1985];
  //sp 112
  sp_rates[111] -= fwd_rates[1985];

  //rxn 1986
  //sp 4
  sp_rates[3] += fwd_rates[1986];
  //sp 6
  sp_rates[5] -= fwd_rates[1986];
  //sp 8
  sp_rates[7] += fwd_rates[1986];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1986];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1986];
  //sp 112
  sp_rates[111] -= fwd_rates[1986];

  //rxn 1987
  //sp 4
  sp_rates[3] += fwd_rates[1987];
  //sp 5
  sp_rates[4] -= fwd_rates[1987];
  //sp 10
  sp_rates[9] += fwd_rates[1987];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1987];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1987];
  //sp 112
  sp_rates[111] -= fwd_rates[1987];

  //rxn 1988
  //sp 4
  sp_rates[3] += fwd_rates[1988];
  //sp 7
  sp_rates[6] += fwd_rates[1988];
  //sp 8
  sp_rates[7] -= fwd_rates[1988];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1988];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1988];
  //sp 112
  sp_rates[111] -= fwd_rates[1988];

  //rxn 1989
  //sp 4
  sp_rates[3] += fwd_rates[1989];
  //sp 9
  sp_rates[8] += fwd_rates[1989];
  //sp 10
  sp_rates[9] -= fwd_rates[1989];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1989];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1989];
  //sp 112
  sp_rates[111] -= fwd_rates[1989];

  //rxn 1990
  //sp 4
  sp_rates[3] += fwd_rates[1990];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1990];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1990];
  //sp 112
  sp_rates[111] -= fwd_rates[1990];
  //sp 22
  sp_rates[21] += fwd_rates[1990];
  //sp 23
  sp_rates[22] -= fwd_rates[1990];

  //rxn 1991
  //sp 4
  sp_rates[3] += fwd_rates[1991];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1991];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1991];
  //sp 112
  sp_rates[111] -= fwd_rates[1991];
  //sp 22
  sp_rates[21] += fwd_rates[1991];
  //sp 24
  sp_rates[23] -= fwd_rates[1991];

  //rxn 1992
  //sp 4
  sp_rates[3] += fwd_rates[1992];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1992];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1992];
  //sp 112
  sp_rates[111] -= fwd_rates[1992];
  //sp 25
  sp_rates[24] += fwd_rates[1992];
  //sp 26
  sp_rates[25] -= fwd_rates[1992];

  //rxn 1993
  //sp 4
  sp_rates[3] += fwd_rates[1993];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1993];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1993];
  //sp 45
  sp_rates[44] += fwd_rates[1993];
  //sp 46
  sp_rates[45] -= fwd_rates[1993];
  //sp 112
  sp_rates[111] -= fwd_rates[1993];

  //rxn 1994
  //sp 4
  sp_rates[3] += fwd_rates[1994];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1994];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1994];
  //sp 45
  sp_rates[44] += fwd_rates[1994];
  //sp 47
  sp_rates[46] -= fwd_rates[1994];
  //sp 112
  sp_rates[111] -= fwd_rates[1994];

  //rxn 1995
  //sp 4
  sp_rates[3] += fwd_rates[1995];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1995];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1995];
  //sp 79
  sp_rates[78] += fwd_rates[1995];
  //sp 112
  sp_rates[111] -= fwd_rates[1995];
  //sp 80
  sp_rates[79] -= fwd_rates[1995];

  //rxn 1996
  //sp 4
  sp_rates[3] += fwd_rates[1996];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1996];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1996];
  //sp 48
  sp_rates[47] += fwd_rates[1996];
  //sp 49
  sp_rates[48] -= fwd_rates[1996];
  //sp 112
  sp_rates[111] -= fwd_rates[1996];

  //rxn 1997
  //sp 4
  sp_rates[3] += fwd_rates[1997];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1997];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1997];
  //sp 14
  sp_rates[13] += fwd_rates[1997];
  //sp 15
  sp_rates[14] -= fwd_rates[1997];
  //sp 112
  sp_rates[111] -= fwd_rates[1997];

  //rxn 1998
  //sp 4
  sp_rates[3] += fwd_rates[1998];
  //sp 107
  sp_rates[106] += 0.85 * fwd_rates[1998];
  //sp 108
  sp_rates[107] += 0.15 * fwd_rates[1998];
  //sp 110
  sp_rates[109] += fwd_rates[1998];
  //sp 111
  sp_rates[110] -= fwd_rates[1998];
  //sp 112
  sp_rates[111] -= fwd_rates[1998];

  //rxn 1999
  //sp 4
  sp_rates[3] += (fwd_rates[1999] - rev_rates[858]) * pres_mod[49];
  //sp 133
  sp_rates[132] = (fwd_rates[1999] - rev_rates[858]) * pres_mod[49];
  //sp 127
  sp_rates[126] = -(fwd_rates[1999] - rev_rates[858]) * pres_mod[49];

  //rxn 2000
  //sp 3
  sp_rates[2] += (fwd_rates[2000] - rev_rates[859]) * pres_mod[50];
  //sp 127
  sp_rates[126] -= (fwd_rates[2000] - rev_rates[859]) * pres_mod[50];
  //sp 131
  sp_rates[130] = (fwd_rates[2000] - rev_rates[859]) * pres_mod[50];

  //rxn 2001
  //sp 3
  sp_rates[2] += (fwd_rates[2001] - rev_rates[860]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2001] - rev_rates[860]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2001] - rev_rates[860]);
  //sp 131
  sp_rates[130] += (fwd_rates[2001] - rev_rates[860]);

  //rxn 2002
  //sp 3
  sp_rates[2] += (fwd_rates[2002] - rev_rates[861]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2002] - rev_rates[861]);
  //sp 133
  sp_rates[132] += (fwd_rates[2002] - rev_rates[861]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2002] - rev_rates[861]);

  //rxn 2003
  //sp 133
  sp_rates[132] += (fwd_rates[2003] - rev_rates[862]);
  //sp 7
  sp_rates[6] += (fwd_rates[2003] - rev_rates[862]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2003] - rev_rates[862]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2003] - rev_rates[862]);

  //rxn 2004
  //sp 133
  sp_rates[132] += (fwd_rates[2004] - rev_rates[863]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2004] - rev_rates[863]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2004] - rev_rates[863]);
  //sp 8
  sp_rates[7] += (fwd_rates[2004] - rev_rates[863]);

  //rxn 2005
  //sp 9
  sp_rates[8] += (fwd_rates[2005] - rev_rates[864]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2005] - rev_rates[864]);
  //sp 133
  sp_rates[132] += (fwd_rates[2005] - rev_rates[864]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2005] - rev_rates[864]);

  //rxn 2006
  //sp 3
  sp_rates[2] += (fwd_rates[2006] - rev_rates[865]);
  //sp 135
  sp_rates[134] = (fwd_rates[2006] - rev_rates[865]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2006] - rev_rates[865]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2006] - rev_rates[865]);

  //rxn 2007
  //sp 117
  sp_rates[116] = (fwd_rates[2007] - rev_rates[866]);
  //sp 4
  sp_rates[3] += (fwd_rates[2007] - rev_rates[866]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2007] - rev_rates[866]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2007] - rev_rates[866]);

  //rxn 2008
  //sp 131
  sp_rates[130] += (fwd_rates[2008] - rev_rates[867]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2008] - rev_rates[867]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2008] - rev_rates[867]);
  //sp 8
  sp_rates[7] += (fwd_rates[2008] - rev_rates[867]);

  //rxn 2009
  //sp 131
  sp_rates[130] += (fwd_rates[2009] - rev_rates[868]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2009] - rev_rates[868]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2009] - rev_rates[868]);
  //sp 8
  sp_rates[7] += (fwd_rates[2009] - rev_rates[868]);

  //rxn 2010
  //sp 117
  sp_rates[116] += (fwd_rates[2010] - rev_rates[869]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2010] - rev_rates[869]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2010] - rev_rates[869]);
  //sp 8
  sp_rates[7] += (fwd_rates[2010] - rev_rates[869]);

  //rxn 2011
  //sp 133
  sp_rates[132] -= (fwd_rates[2011] - rev_rates[870]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2011] - rev_rates[870]);
  //sp 6
  sp_rates[5] += (fwd_rates[2011] - rev_rates[870]);
  //sp 134
  sp_rates[133] = (fwd_rates[2011] - rev_rates[870]);

  //rxn 2012
  //sp 10
  sp_rates[9] -= (fwd_rates[2012] - rev_rates[871]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2012] - rev_rates[871]);
  //sp 134
  sp_rates[133] += (fwd_rates[2012] - rev_rates[871]);
  //sp 8
  sp_rates[7] += (fwd_rates[2012] - rev_rates[871]);

  //rxn 2013
  //sp 10
  sp_rates[9] -= (fwd_rates[2013] - rev_rates[872]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2013] - rev_rates[872]);
  //sp 5
  sp_rates[4] += (fwd_rates[2013] - rev_rates[872]);
  //sp 127
  sp_rates[126] += (fwd_rates[2013] - rev_rates[872]);

  //rxn 2014
  //sp 121
  sp_rates[120] = (fwd_rates[2014] - rev_rates[873]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2014] - rev_rates[873]);
  //sp 4
  sp_rates[3] += (fwd_rates[2014] - rev_rates[873]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2014] - rev_rates[873]);

  //rxn 2015
  //sp 129
  sp_rates[128] = (fwd_rates[2015] - rev_rates[874]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2015] - rev_rates[874]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2015] - rev_rates[874]);
  //sp 127
  sp_rates[126] += (fwd_rates[2015] - rev_rates[874]);

  //rxn 2016
  //sp 129
  sp_rates[128] += (fwd_rates[2016] - rev_rates[875]);
  //sp 131
  sp_rates[130] -= 2.0 * (fwd_rates[2016] - rev_rates[875]);
  //sp 133
  sp_rates[132] += (fwd_rates[2016] - rev_rates[875]);

  //rxn 2017
  //sp 121
  sp_rates[120] += (fwd_rates[2017] - rev_rates[876]);
  //sp 131
  sp_rates[130] -= 2.0 * (fwd_rates[2017] - rev_rates[876]);

  //rxn 2018
  //sp 121
  sp_rates[120] += (fwd_rates[2018] - rev_rates[877]);
  //sp 3
  sp_rates[2] += (fwd_rates[2018] - rev_rates[877]);
  //sp 133
  sp_rates[132] -= 2.0 * (fwd_rates[2018] - rev_rates[877]);

  //rxn 2019
  //sp 131
  sp_rates[130] += (fwd_rates[2019] - rev_rates[878]);
  //sp 133
  sp_rates[132] -= 2.0 * (fwd_rates[2019] - rev_rates[878]);
  //sp 127
  sp_rates[126] += (fwd_rates[2019] - rev_rates[878]);

  //rxn 2020
  //sp 122
  sp_rates[121] = (fwd_rates[2020] - rev_rates[879]);
  //sp 3
  sp_rates[2] += (fwd_rates[2020] - rev_rates[879]);
  //sp 133
  sp_rates[132] -= 2.0 * (fwd_rates[2020] - rev_rates[879]);

  //rxn 2021
  //sp 133
  sp_rates[132] -= 2.0 * (fwd_rates[2021] - rev_rates[880]);
  //sp 128
  sp_rates[127] = (fwd_rates[2021] - rev_rates[880]);

  //rxn 2022
  //sp 4
  sp_rates[3] += (fwd_rates[2022] - rev_rates[881]);
  //sp 133
  sp_rates[132] -= 2.0 * (fwd_rates[2022] - rev_rates[881]);
  //sp 135
  sp_rates[134] += (fwd_rates[2022] - rev_rates[881]);

  //rxn 2023
  //sp 133
  sp_rates[132] += (fwd_rates[2023] - rev_rates[882]) * pres_mod[51];
  //sp 125
  sp_rates[124] = -(fwd_rates[2023] - rev_rates[882]) * pres_mod[51];
  //sp 8
  sp_rates[7] += (fwd_rates[2023] - rev_rates[882]) * pres_mod[51];

  //rxn 2024
  //sp 3
  sp_rates[2] += (fwd_rates[2024] - rev_rates[883]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2024] - rev_rates[883]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2024] - rev_rates[883]);
  //sp 126
  sp_rates[125] = (fwd_rates[2024] - rev_rates[883]);

  //rxn 2025
  //sp 3
  sp_rates[2] += (fwd_rates[2025] - rev_rates[884]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2025] - rev_rates[884]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2025] - rev_rates[884]);
  //sp 134
  sp_rates[133] += (fwd_rates[2025] - rev_rates[884]);

  //rxn 2026
  //sp 126
  sp_rates[125] += (fwd_rates[2026] - rev_rates[885]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2026] - rev_rates[885]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2026] - rev_rates[885]);
  //sp 8
  sp_rates[7] += (fwd_rates[2026] - rev_rates[885]);

  //rxn 2027
  //sp 125
  sp_rates[124] -= (fwd_rates[2027] - rev_rates[886]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2027] - rev_rates[886]);
  //sp 134
  sp_rates[133] += (fwd_rates[2027] - rev_rates[886]);
  //sp 8
  sp_rates[7] += (fwd_rates[2027] - rev_rates[886]);

  //rxn 2028
  //sp 125
  sp_rates[124] -= (fwd_rates[2028] - rev_rates[887]);
  //sp 126
  sp_rates[125] += (fwd_rates[2028] - rev_rates[887]);
  //sp 7
  sp_rates[6] += (fwd_rates[2028] - rev_rates[887]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2028] - rev_rates[887]);

  //rxn 2029
  //sp 125
  sp_rates[124] -= (fwd_rates[2029] - rev_rates[888]);
  //sp 134
  sp_rates[133] += (fwd_rates[2029] - rev_rates[888]);
  //sp 7
  sp_rates[6] += (fwd_rates[2029] - rev_rates[888]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2029] - rev_rates[888]);

  //rxn 2030
  //sp 133
  sp_rates[132] -= (fwd_rates[2030] - rev_rates[889]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2030] - rev_rates[889]);
  //sp 126
  sp_rates[125] += (fwd_rates[2030] - rev_rates[889]);
  //sp 127
  sp_rates[126] += (fwd_rates[2030] - rev_rates[889]);

  //rxn 2031
  //sp 133
  sp_rates[132] -= (fwd_rates[2031] - rev_rates[890]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2031] - rev_rates[890]);
  //sp 134
  sp_rates[133] += (fwd_rates[2031] - rev_rates[890]);
  //sp 127
  sp_rates[126] += (fwd_rates[2031] - rev_rates[890]);

  //rxn 2032
  //sp 133
  sp_rates[132] += (fwd_rates[2032] - rev_rates[891]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2032] - rev_rates[891]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2032] - rev_rates[891]);
  //sp 126
  sp_rates[125] += (fwd_rates[2032] - rev_rates[891]);

  //rxn 2033
  //sp 133
  sp_rates[132] += (fwd_rates[2033] - rev_rates[892]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2033] - rev_rates[892]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2033] - rev_rates[892]);
  //sp 134
  sp_rates[133] += (fwd_rates[2033] - rev_rates[892]);

  //rxn 2034
  //sp 9
  sp_rates[8] += (fwd_rates[2034] - rev_rates[893]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2034] - rev_rates[893]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2034] - rev_rates[893]);
  //sp 126
  sp_rates[125] += (fwd_rates[2034] - rev_rates[893]);

  //rxn 2035
  //sp 9
  sp_rates[8] += (fwd_rates[2035] - rev_rates[894]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2035] - rev_rates[894]);
  //sp 125
  sp_rates[124] -= (fwd_rates[2035] - rev_rates[894]);
  //sp 134
  sp_rates[133] += (fwd_rates[2035] - rev_rates[894]);

  //rxn 2036
  //sp 4
  sp_rates[3] += (fwd_rates[2036] - rev_rates[895]) * pres_mod[52];
  //sp 117
  sp_rates[116] += (fwd_rates[2036] - rev_rates[895]) * pres_mod[52];
  //sp 126
  sp_rates[125] -= (fwd_rates[2036] - rev_rates[895]) * pres_mod[52];

  //rxn 2037
  //sp 117
  sp_rates[116] += (fwd_rates[2037] - rev_rates[896]);
  //sp 10
  sp_rates[9] += (fwd_rates[2037] - rev_rates[896]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2037] - rev_rates[896]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2037] - rev_rates[896]);

  //rxn 2038
  //sp 4
  sp_rates[3] -= (fwd_rates[2038] - rev_rates[897]);
  //sp 133
  sp_rates[132] += (fwd_rates[2038] - rev_rates[897]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2038] - rev_rates[897]);
  //sp 8
  sp_rates[7] += (fwd_rates[2038] - rev_rates[897]);

  //rxn 2039
  //sp 3
  sp_rates[2] += (fwd_rates[2039] - rev_rates[898]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2039] - rev_rates[898]);
  //sp 117
  sp_rates[116] += (fwd_rates[2039] - rev_rates[898]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2039] - rev_rates[898]);

  //rxn 2040
  //sp 126
  sp_rates[125] -= (fwd_rates[2040] - rev_rates[899]);
  //sp 117
  sp_rates[116] += (fwd_rates[2040] - rev_rates[899]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2040] - rev_rates[899]);
  //sp 8
  sp_rates[7] += (fwd_rates[2040] - rev_rates[899]);

  //rxn 2041
  //sp 126
  sp_rates[125] -= (fwd_rates[2041] - rev_rates[900]);
  //sp 117
  sp_rates[116] += (fwd_rates[2041] - rev_rates[900]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2041] - rev_rates[900]);
  //sp 8
  sp_rates[7] += (fwd_rates[2041] - rev_rates[900]);

  //rxn 2042
  //sp 117
  sp_rates[116] += (fwd_rates[2042] - rev_rates[901]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2042] - rev_rates[901]);
  //sp 7
  sp_rates[6] += (fwd_rates[2042] - rev_rates[901]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2042] - rev_rates[901]);

  //rxn 2043
  //sp 9
  sp_rates[8] += (fwd_rates[2043] - rev_rates[902]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2043] - rev_rates[902]);
  //sp 117
  sp_rates[116] += (fwd_rates[2043] - rev_rates[902]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2043] - rev_rates[902]);

  //rxn 2044
  //sp 125
  sp_rates[124] += (fwd_rates[2044] - rev_rates[903]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2044] - rev_rates[903]);
  //sp 5
  sp_rates[4] += (fwd_rates[2044] - rev_rates[903]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2044] - rev_rates[903]);

  //rxn 2045
  //sp 117
  sp_rates[116] += (fwd_rates[2045] - rev_rates[904]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2045] - rev_rates[904]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2045] - rev_rates[904]);
  //sp 127
  sp_rates[126] += (fwd_rates[2045] - rev_rates[904]);

  //rxn 2046
  //sp 133
  sp_rates[132] -= (fwd_rates[2046] - rev_rates[905]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2046] - rev_rates[905]);
  //sp 135
  sp_rates[134] += (fwd_rates[2046] - rev_rates[905]);
  //sp 8
  sp_rates[7] += (fwd_rates[2046] - rev_rates[905]);

  //rxn 2047
  //sp 122
  sp_rates[121] += (fwd_rates[2047] - rev_rates[906]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2047] - rev_rates[906]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2047] - rev_rates[906]);
  //sp 7
  sp_rates[6] += (fwd_rates[2047] - rev_rates[906]);

  //rxn 2048
  //sp 116
  sp_rates[115] = -(fwd_rates[2048] - rev_rates[907]);
  //sp 117
  sp_rates[116] += (fwd_rates[2048] - rev_rates[907]);
  //sp 126
  sp_rates[125] -= (fwd_rates[2048] - rev_rates[907]);
  //sp 119
  sp_rates[118] = (fwd_rates[2048] - rev_rates[907]);

  //rxn 2049
  //sp 115
  sp_rates[114] = (fwd_rates[2049] - rev_rates[908]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2049] - rev_rates[908]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2049] - rev_rates[908]);
  //sp 7
  sp_rates[6] += (fwd_rates[2049] - rev_rates[908]);

  //rxn 2050
  //sp 114
  sp_rates[113] = (fwd_rates[2050] - rev_rates[909]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2050] - rev_rates[909]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2050] - rev_rates[909]);
  //sp 134
  sp_rates[133] += (fwd_rates[2050] - rev_rates[909]);

  //rxn 2051
  //sp 114
  sp_rates[113] -= (fwd_rates[2051] - rev_rates[910]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2051] - rev_rates[910]);
  //sp 1
  (*dy_N) = (fwd_rates[2051] - rev_rates[910]);
  //sp 7
  sp_rates[6] += (fwd_rates[2051] - rev_rates[910]);

  //rxn 2052
  //sp 114
  sp_rates[113] -= (fwd_rates[2052] - rev_rates[911]);
  //sp 132
  sp_rates[131] = (fwd_rates[2052] - rev_rates[911]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2052] - rev_rates[911]);
  //sp 8
  sp_rates[7] += (fwd_rates[2052] - rev_rates[911]);

  //rxn 2053
  //sp 114
  sp_rates[113] -= (fwd_rates[2053] - rev_rates[912]);
  //sp 124
  sp_rates[123] = (fwd_rates[2053] - rev_rates[912]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2053] - rev_rates[912]);

  //rxn 2054
  //sp 3
  sp_rates[2] += (fwd_rates[2054] - rev_rates[913]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2054] - rev_rates[913]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2054] - rev_rates[913]);
  //sp 129
  sp_rates[128] += (fwd_rates[2054] - rev_rates[913]);

  //rxn 2055
  //sp 114
  sp_rates[113] += (fwd_rates[2055] - rev_rates[914]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2055] - rev_rates[914]);
  //sp 4
  sp_rates[3] += (fwd_rates[2055] - rev_rates[914]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2055] - rev_rates[914]);

  //rxn 2056
  //sp 131
  sp_rates[130] -= (fwd_rates[2056] - rev_rates[915]);
  //sp 4
  sp_rates[3] += (fwd_rates[2056] - rev_rates[915]);
  //sp 117
  sp_rates[116] += (fwd_rates[2056] - rev_rates[915]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2056] - rev_rates[915]);

  //rxn 2057
  //sp 129
  sp_rates[128] += (fwd_rates[2057] - rev_rates[916]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2057] - rev_rates[916]);
  //sp 7
  sp_rates[6] += (fwd_rates[2057] - rev_rates[916]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2057] - rev_rates[916]);

  //rxn 2058
  //sp 114
  sp_rates[113] += (fwd_rates[2058] - rev_rates[917]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2058] - rev_rates[917]);
  //sp 3
  sp_rates[2] += (fwd_rates[2058] - rev_rates[917]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2058] - rev_rates[917]);

  //rxn 2059
  //sp 131
  sp_rates[130] += (fwd_rates[2059] - rev_rates[918]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2059] - rev_rates[918]);
  //sp 7
  sp_rates[6] += (fwd_rates[2059] - rev_rates[918]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2059] - rev_rates[918]);

  //rxn 2060
  //sp 117
  sp_rates[116] += (fwd_rates[2060] - rev_rates[919]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2060] - rev_rates[919]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2060] - rev_rates[919]);
  //sp 6
  sp_rates[5] += (fwd_rates[2060] - rev_rates[919]);

  //rxn 2061
  //sp 114
  sp_rates[113] += (fwd_rates[2061] - rev_rates[920]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2061] - rev_rates[920]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2061] - rev_rates[920]);
  //sp 8
  sp_rates[7] += (fwd_rates[2061] - rev_rates[920]);

  //rxn 2062
  //sp 131
  sp_rates[130] -= 2.0 * (fwd_rates[2062] - rev_rates[921]);
  //sp 4
  sp_rates[3] += 2.0 * (fwd_rates[2062] - rev_rates[921]);
  //sp 1
  (*dy_N) += (fwd_rates[2062] - rev_rates[921]);

  //rxn 2063
  //sp 129
  sp_rates[128] -= (fwd_rates[2063] - rev_rates[922]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2063] - rev_rates[922]);
  //sp 4
  sp_rates[3] += (fwd_rates[2063] - rev_rates[922]);
  //sp 1
  (*dy_N) += (fwd_rates[2063] - rev_rates[922]);

  //rxn 2064
  //sp 114
  sp_rates[113] -= (fwd_rates[2064] - rev_rates[923]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2064] - rev_rates[923]);
  //sp 1
  (*dy_N) += (fwd_rates[2064] - rev_rates[923]);
  //sp 8
  sp_rates[7] += (fwd_rates[2064] - rev_rates[923]);

  //rxn 2065
  //sp 114
  sp_rates[113] -= (fwd_rates[2065] - rev_rates[924]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2065] - rev_rates[924]);
  //sp 4
  sp_rates[3] += (fwd_rates[2065] - rev_rates[924]);
  //sp 115
  sp_rates[114] += (fwd_rates[2065] - rev_rates[924]);

  //rxn 2066
  //sp 131
  sp_rates[130] -= (fwd_rates[2066] - rev_rates[925]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2066] - rev_rates[925]);
  //sp 115
  sp_rates[114] += (fwd_rates[2066] - rev_rates[925]);
  //sp 8
  sp_rates[7] += (fwd_rates[2066] - rev_rates[925]);

  //rxn 2067
  //sp 114
  sp_rates[113] += (fwd_rates[2067] - rev_rates[926]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2067] - rev_rates[926]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2067] - rev_rates[926]);
  //sp 117
  sp_rates[116] += (fwd_rates[2067] - rev_rates[926]);

  //rxn 2068
  //sp 129
  sp_rates[128] -= (fwd_rates[2068] - rev_rates[927]);
  //sp 114
  sp_rates[113] += (fwd_rates[2068] - rev_rates[927]);
  //sp 4
  sp_rates[3] += (fwd_rates[2068] - rev_rates[927]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2068] - rev_rates[927]);

  //rxn 2069
  //sp 129
  sp_rates[128] -= (fwd_rates[2069] - rev_rates[928]);
  //sp 114
  sp_rates[113] += (fwd_rates[2069] - rev_rates[928]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2069] - rev_rates[928]);
  //sp 6
  sp_rates[5] += (fwd_rates[2069] - rev_rates[928]);

  //rxn 2070
  //sp 129
  sp_rates[128] -= (fwd_rates[2070] - rev_rates[929]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2070] - rev_rates[929]);
  //sp 6
  sp_rates[5] += (fwd_rates[2070] - rev_rates[929]);
  //sp 1
  (*dy_N) += (fwd_rates[2070] - rev_rates[929]);

  //rxn 2071
  //sp 122
  sp_rates[121] += (fwd_rates[2071] - rev_rates[930]);
  //sp 3
  sp_rates[2] += (fwd_rates[2071] - rev_rates[930]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2071] - rev_rates[930]);

  //rxn 2072
  //sp 3
  sp_rates[2] += (fwd_rates[2072] - rev_rates[931]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2072] - rev_rates[931]);
  //sp 135
  sp_rates[134] += (fwd_rates[2072] - rev_rates[931]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2072] - rev_rates[931]);

  //rxn 2073
  //sp 8
  sp_rates[7] += (fwd_rates[2073] - rev_rates[932]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2073] - rev_rates[932]);
  //sp 135
  sp_rates[134] += (fwd_rates[2073] - rev_rates[932]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2073] - rev_rates[932]);

  //rxn 2074
  //sp 121
  sp_rates[120] += (fwd_rates[2074] - rev_rates[933]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2074] - rev_rates[933]);
  //sp 7
  sp_rates[6] += (fwd_rates[2074] - rev_rates[933]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2074] - rev_rates[933]);

  //rxn 2075
  //sp 135
  sp_rates[134] += (fwd_rates[2075] - rev_rates[934]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2075] - rev_rates[934]);
  //sp 7
  sp_rates[6] += (fwd_rates[2075] - rev_rates[934]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2075] - rev_rates[934]);

  //rxn 2076
  //sp 9
  sp_rates[8] += (fwd_rates[2076] - rev_rates[935]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2076] - rev_rates[935]);
  //sp 135
  sp_rates[134] += (fwd_rates[2076] - rev_rates[935]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2076] - rev_rates[935]);

  //rxn 2077
  //sp 135
  sp_rates[134] += (fwd_rates[2077] - rev_rates[936]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2077] - rev_rates[936]);
  //sp 127
  sp_rates[126] += (fwd_rates[2077] - rev_rates[936]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2077] - rev_rates[936]);

  //rxn 2078
  //sp 135
  sp_rates[134] += (fwd_rates[2078] - rev_rates[937]);
  //sp 14
  sp_rates[13] += (fwd_rates[2078] - rev_rates[937]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2078] - rev_rates[937]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2078] - rev_rates[937]);

  //rxn 2079
  //sp 114
  sp_rates[113] -= (fwd_rates[2079] - rev_rates[938]);
  //sp 117
  sp_rates[116] += (fwd_rates[2079] - rev_rates[938]);
  //sp 135
  sp_rates[134] += (fwd_rates[2079] - rev_rates[938]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2079] - rev_rates[938]);

  //rxn 2080
  //sp 135
  sp_rates[134] += (fwd_rates[2080] - rev_rates[939]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2080] - rev_rates[939]);
  //sp 119
  sp_rates[118] += (fwd_rates[2080] - rev_rates[939]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2080] - rev_rates[939]);

  //rxn 2081
  //sp 116
  sp_rates[115] -= (fwd_rates[2081] - rev_rates[940]);
  //sp 118
  sp_rates[117] = (fwd_rates[2081] - rev_rates[940]);
  //sp 135
  sp_rates[134] += (fwd_rates[2081] - rev_rates[940]);
  //sp 128
  sp_rates[127] -= (fwd_rates[2081] - rev_rates[940]);

  //rxn 2082
  //sp 121
  sp_rates[120] += (fwd_rates[2082] - rev_rates[941]);
  //sp 4
  sp_rates[3] += (fwd_rates[2082] - rev_rates[941]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2082] - rev_rates[941]);

  //rxn 2083
  //sp 121
  sp_rates[120] += (fwd_rates[2083] - rev_rates[942]);
  //sp 3
  sp_rates[2] += (fwd_rates[2083] - rev_rates[942]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2083] - rev_rates[942]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2083] - rev_rates[942]);

  //rxn 2084
  //sp 133
  sp_rates[132] += (fwd_rates[2084] - rev_rates[943]);
  //sp 117
  sp_rates[116] += (fwd_rates[2084] - rev_rates[943]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2084] - rev_rates[943]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2084] - rev_rates[943]);

  //rxn 2085
  //sp 121
  sp_rates[120] += (fwd_rates[2085] - rev_rates[944]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2085] - rev_rates[944]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2085] - rev_rates[944]);
  //sp 8
  sp_rates[7] += (fwd_rates[2085] - rev_rates[944]);

  //rxn 2086
  //sp 124
  sp_rates[123] += (fwd_rates[2086] - rev_rates[945]);
  //sp 4
  sp_rates[3] += (fwd_rates[2086] - rev_rates[945]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2086] - rev_rates[945]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2086] - rev_rates[945]);

  //rxn 2087
  //sp 121
  sp_rates[120] += (fwd_rates[2087] - rev_rates[946]);
  //sp 7
  sp_rates[6] += (fwd_rates[2087] - rev_rates[946]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2087] - rev_rates[946]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2087] - rev_rates[946]);

  //rxn 2088
  //sp 122
  sp_rates[121] += (fwd_rates[2088] - rev_rates[947]);
  //sp 7
  sp_rates[6] += (fwd_rates[2088] - rev_rates[947]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2088] - rev_rates[947]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2088] - rev_rates[947]);

  //rxn 2089
  //sp 122
  sp_rates[121] += (fwd_rates[2089] - rev_rates[948]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2089] - rev_rates[948]);
  //sp 14
  sp_rates[13] += (fwd_rates[2089] - rev_rates[948]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2089] - rev_rates[948]);

  //rxn 2090
  //sp 121
  sp_rates[120] += (fwd_rates[2090] - rev_rates[949]);
  //sp 127
  sp_rates[126] += (fwd_rates[2090] - rev_rates[949]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2090] - rev_rates[949]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2090] - rev_rates[949]);

  //rxn 2091
  //sp 122
  sp_rates[121] += (fwd_rates[2091] - rev_rates[950]);
  //sp 127
  sp_rates[126] += (fwd_rates[2091] - rev_rates[950]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2091] - rev_rates[950]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2091] - rev_rates[950]);

  //rxn 2092
  //sp 9
  sp_rates[8] += (fwd_rates[2092] - rev_rates[951]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2092] - rev_rates[951]);
  //sp 121
  sp_rates[120] += (fwd_rates[2092] - rev_rates[951]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2092] - rev_rates[951]);

  //rxn 2093
  //sp 10
  sp_rates[9] -= (fwd_rates[2093] - rev_rates[952]);
  //sp 5
  sp_rates[4] += (fwd_rates[2093] - rev_rates[952]);
  //sp 135
  sp_rates[134] -= (fwd_rates[2093] - rev_rates[952]);
  //sp 128
  sp_rates[127] += (fwd_rates[2093] - rev_rates[952]);

  //rxn 2094
  //sp 121
  sp_rates[120] -= (fwd_rates[2094] - rev_rates[953]);
  //sp 4
  sp_rates[3] += (fwd_rates[2094] - rev_rates[953]);
  //sp 132
  sp_rates[131] += (fwd_rates[2094] - rev_rates[953]);

  //rxn 2095
  //sp 121
  sp_rates[120] -= (fwd_rates[2095] - rev_rates[954]);
  //sp 4
  sp_rates[3] += (fwd_rates[2095] - rev_rates[954]);
  //sp 132
  sp_rates[131] += (fwd_rates[2095] - rev_rates[954]);

  //rxn 2096
  //sp 121
  sp_rates[120] -= (fwd_rates[2096] - rev_rates[955]);
  //sp 3
  sp_rates[2] += (fwd_rates[2096] - rev_rates[955]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2096] - rev_rates[955]);
  //sp 132
  sp_rates[131] += (fwd_rates[2096] - rev_rates[955]);

  //rxn 2097
  //sp 121
  sp_rates[120] -= (fwd_rates[2097] - rev_rates[956]);
  //sp 132
  sp_rates[131] += (fwd_rates[2097] - rev_rates[956]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2097] - rev_rates[956]);
  //sp 8
  sp_rates[7] += (fwd_rates[2097] - rev_rates[956]);

  //rxn 2098
  //sp 121
  sp_rates[120] -= (fwd_rates[2098] - rev_rates[957]);
  //sp 132
  sp_rates[131] += (fwd_rates[2098] - rev_rates[957]);
  //sp 7
  sp_rates[6] += (fwd_rates[2098] - rev_rates[957]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2098] - rev_rates[957]);

  //rxn 2099
  //sp 121
  sp_rates[120] -= (fwd_rates[2099] - rev_rates[958]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2099] - rev_rates[958]);
  //sp 115
  sp_rates[114] += (fwd_rates[2099] - rev_rates[958]);
  //sp 133
  sp_rates[132] += (fwd_rates[2099] - rev_rates[958]);

  //rxn 2100
  //sp 121
  sp_rates[120] -= (fwd_rates[2100] - rev_rates[959]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2100] - rev_rates[959]);
  //sp 132
  sp_rates[131] += (fwd_rates[2100] - rev_rates[959]);
  //sp 9
  sp_rates[8] += (fwd_rates[2100] - rev_rates[959]);

  //rxn 2101
  //sp 121
  sp_rates[120] -= (fwd_rates[2101] - rev_rates[960]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2101] - rev_rates[960]);
  //sp 132
  sp_rates[131] += (fwd_rates[2101] - rev_rates[960]);
  //sp 133
  sp_rates[132] += (fwd_rates[2101] - rev_rates[960]);

  //rxn 2102
  //sp 121
  sp_rates[120] -= (fwd_rates[2102] - rev_rates[961]);
  //sp 132
  sp_rates[131] += (fwd_rates[2102] - rev_rates[961]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2102] - rev_rates[961]);
  //sp 127
  sp_rates[126] += (fwd_rates[2102] - rev_rates[961]);

  //rxn 2103
  //sp 121
  sp_rates[120] -= (fwd_rates[2103] - rev_rates[962]);
  //sp 122
  sp_rates[121] += (fwd_rates[2103] - rev_rates[962]);

  //rxn 2104
  //sp 122
  sp_rates[121] -= (fwd_rates[2104] - rev_rates[963]);
  //sp 4
  sp_rates[3] += (fwd_rates[2104] - rev_rates[963]);
  //sp 132
  sp_rates[131] += (fwd_rates[2104] - rev_rates[963]);

  //rxn 2105
  //sp 122
  sp_rates[121] -= (fwd_rates[2105] - rev_rates[964]);
  //sp 4
  sp_rates[3] += (fwd_rates[2105] - rev_rates[964]);
  //sp 132
  sp_rates[131] += (fwd_rates[2105] - rev_rates[964]);

  //rxn 2106
  //sp 122
  sp_rates[121] -= (fwd_rates[2106] - rev_rates[965]);
  //sp 116
  sp_rates[115] += (fwd_rates[2106] - rev_rates[965]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2106] - rev_rates[965]);
  //sp 133
  sp_rates[132] += (fwd_rates[2106] - rev_rates[965]);

  //rxn 2107
  //sp 121
  sp_rates[120] += (fwd_rates[2107] - rev_rates[966]);
  //sp 122
  sp_rates[121] -= (fwd_rates[2107] - rev_rates[966]);

  //rxn 2108
  //sp 122
  sp_rates[121] -= (fwd_rates[2108] - rev_rates[967]);
  //sp 3
  sp_rates[2] += (fwd_rates[2108] - rev_rates[967]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2108] - rev_rates[967]);
  //sp 132
  sp_rates[131] += (fwd_rates[2108] - rev_rates[967]);

  //rxn 2109
  //sp 122
  sp_rates[121] -= (fwd_rates[2109] - rev_rates[968]);
  //sp 133
  sp_rates[132] += (fwd_rates[2109] - rev_rates[968]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2109] - rev_rates[968]);
  //sp 114
  sp_rates[113] += (fwd_rates[2109] - rev_rates[968]);

  //rxn 2110
  //sp 122
  sp_rates[121] -= (fwd_rates[2110] - rev_rates[969]);
  //sp 132
  sp_rates[131] += (fwd_rates[2110] - rev_rates[969]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2110] - rev_rates[969]);
  //sp 8
  sp_rates[7] += (fwd_rates[2110] - rev_rates[969]);

  //rxn 2111
  //sp 124
  sp_rates[123] += (fwd_rates[2111] - rev_rates[970]);
  //sp 122
  sp_rates[121] -= (fwd_rates[2111] - rev_rates[970]);
  //sp 4
  sp_rates[3] += (fwd_rates[2111] - rev_rates[970]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2111] - rev_rates[970]);

  //rxn 2112
  //sp 122
  sp_rates[121] -= (fwd_rates[2112] - rev_rates[971]);
  //sp 132
  sp_rates[131] += (fwd_rates[2112] - rev_rates[971]);
  //sp 7
  sp_rates[6] += (fwd_rates[2112] - rev_rates[971]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2112] - rev_rates[971]);

  //rxn 2113
  //sp 122
  sp_rates[121] -= (fwd_rates[2113] - rev_rates[972]);
  //sp 132
  sp_rates[131] += (fwd_rates[2113] - rev_rates[972]);
  //sp 14
  sp_rates[13] += (fwd_rates[2113] - rev_rates[972]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2113] - rev_rates[972]);

  //rxn 2114
  //sp 122
  sp_rates[121] -= (fwd_rates[2114] - rev_rates[973]);
  //sp 132
  sp_rates[131] += (fwd_rates[2114] - rev_rates[973]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2114] - rev_rates[973]);
  //sp 127
  sp_rates[126] += (fwd_rates[2114] - rev_rates[973]);

  //rxn 2115
  //sp 10
  sp_rates[9] -= (fwd_rates[2115] - rev_rates[974]);
  //sp 124
  sp_rates[123] += (fwd_rates[2115] - rev_rates[974]);
  //sp 122
  sp_rates[121] -= (fwd_rates[2115] - rev_rates[974]);
  //sp 8
  sp_rates[7] += (fwd_rates[2115] - rev_rates[974]);

  //rxn 2116
  //sp 9
  sp_rates[8] += (fwd_rates[2116] - rev_rates[975]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2116] - rev_rates[975]);
  //sp 132
  sp_rates[131] += (fwd_rates[2116] - rev_rates[975]);
  //sp 122
  sp_rates[121] -= (fwd_rates[2116] - rev_rates[975]);

  //rxn 2117
  //sp 124
  sp_rates[123] -= (fwd_rates[2117] - rev_rates[976]);
  //sp 1
  (*dy_N) += (fwd_rates[2117] - rev_rates[976]);
  //sp 7
  sp_rates[6] += (fwd_rates[2117] - rev_rates[976]);

  //rxn 2118
  //sp 115
  sp_rates[114] -= (fwd_rates[2118] - rev_rates[977]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2118] - rev_rates[977]);
  //sp 123
  sp_rates[122] = (fwd_rates[2118] - rev_rates[977]);

  //rxn 2119
  //sp 124
  sp_rates[123] -= (fwd_rates[2119] - rev_rates[978]);
  //sp 3
  sp_rates[2] += (fwd_rates[2119] - rev_rates[978]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2119] - rev_rates[978]);
  //sp 123
  sp_rates[122] += (fwd_rates[2119] - rev_rates[978]);

  //rxn 2120
  //sp 123
  sp_rates[122] += (fwd_rates[2120] - rev_rates[979]);
  //sp 124
  sp_rates[123] -= (fwd_rates[2120] - rev_rates[979]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2120] - rev_rates[979]);
  //sp 8
  sp_rates[7] += (fwd_rates[2120] - rev_rates[979]);

  //rxn 2121
  //sp 123
  sp_rates[122] += (fwd_rates[2121] - rev_rates[980]);
  //sp 124
  sp_rates[123] -= (fwd_rates[2121] - rev_rates[980]);
  //sp 7
  sp_rates[6] += (fwd_rates[2121] - rev_rates[980]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2121] - rev_rates[980]);

  //rxn 2122
  //sp 123
  sp_rates[122] += (fwd_rates[2122] - rev_rates[981]);
  //sp 124
  sp_rates[123] -= (fwd_rates[2122] - rev_rates[981]);
  //sp 14
  sp_rates[13] += (fwd_rates[2122] - rev_rates[981]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2122] - rev_rates[981]);

  //rxn 2123
  //sp 123
  sp_rates[122] += (fwd_rates[2123] - rev_rates[982]);
  //sp 124
  sp_rates[123] -= (fwd_rates[2123] - rev_rates[982]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2123] - rev_rates[982]);
  //sp 127
  sp_rates[126] += (fwd_rates[2123] - rev_rates[982]);

  //rxn 2124
  //sp 9
  sp_rates[8] += (fwd_rates[2124] - rev_rates[983]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2124] - rev_rates[983]);
  //sp 123
  sp_rates[122] += (fwd_rates[2124] - rev_rates[983]);
  //sp 124
  sp_rates[123] -= (fwd_rates[2124] - rev_rates[983]);

  //rxn 2125
  //sp 4
  sp_rates[3] += (fwd_rates[2125] - rev_rates[984]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2125] - rev_rates[984]);
  //sp 1
  (*dy_N) += (fwd_rates[2125] - rev_rates[984]);

  //rxn 2126
  //sp 3
  sp_rates[2] += (fwd_rates[2126] - rev_rates[985]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2126] - rev_rates[985]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2126] - rev_rates[985]);
  //sp 1
  (*dy_N) += (fwd_rates[2126] - rev_rates[985]);

  //rxn 2127
  //sp 4
  sp_rates[3] += (fwd_rates[2127] - rev_rates[986]);
  //sp 115
  sp_rates[114] += (fwd_rates[2127] - rev_rates[986]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2127] - rev_rates[986]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2127] - rev_rates[986]);

  //rxn 2128
  //sp 114
  sp_rates[113] += (fwd_rates[2128] - rev_rates[987]);
  //sp 131
  sp_rates[130] += (fwd_rates[2128] - rev_rates[987]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2128] - rev_rates[987]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2128] - rev_rates[987]);

  //rxn 2129
  //sp 132
  sp_rates[131] -= (fwd_rates[2129] - rev_rates[988]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2129] - rev_rates[988]);
  //sp 1
  (*dy_N) += (fwd_rates[2129] - rev_rates[988]);
  //sp 8
  sp_rates[7] += (fwd_rates[2129] - rev_rates[988]);

  //rxn 2130
  //sp 132
  sp_rates[131] -= (fwd_rates[2130] - rev_rates[989]);
  //sp 7
  sp_rates[6] += (fwd_rates[2130] - rev_rates[989]);
  //sp 1
  (*dy_N) += (fwd_rates[2130] - rev_rates[989]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2130] - rev_rates[989]);

  //rxn 2131
  //sp 10
  sp_rates[9] += (fwd_rates[2131] - rev_rates[990]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2131] - rev_rates[990]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2131] - rev_rates[990]);
  //sp 1
  (*dy_N) += (fwd_rates[2131] - rev_rates[990]);

  //rxn 2132
  //sp 131
  sp_rates[130] -= (fwd_rates[2132] - rev_rates[991]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2132] - rev_rates[991]);
  //sp 133
  sp_rates[132] += (fwd_rates[2132] - rev_rates[991]);
  //sp 1
  (*dy_N) += (fwd_rates[2132] - rev_rates[991]);

  //rxn 2133
  //sp 127
  sp_rates[126] += (fwd_rates[2133] - rev_rates[992]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2133] - rev_rates[992]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2133] - rev_rates[992]);
  //sp 1
  (*dy_N) += (fwd_rates[2133] - rev_rates[992]);

  //rxn 2134
  //sp 9
  sp_rates[8] += (fwd_rates[2134] - rev_rates[993]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2134] - rev_rates[993]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2134] - rev_rates[993]);
  //sp 1
  (*dy_N) += (fwd_rates[2134] - rev_rates[993]);

  //rxn 2135
  //sp 10
  sp_rates[9] -= (fwd_rates[2135] - rev_rates[994]);
  //sp 123
  sp_rates[122] += (fwd_rates[2135] - rev_rates[994]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2135] - rev_rates[994]);
  //sp 8
  sp_rates[7] += (fwd_rates[2135] - rev_rates[994]);

  //rxn 2136
  //sp 114
  sp_rates[113] -= (fwd_rates[2136] - rev_rates[995]);
  //sp 132
  sp_rates[131] -= (fwd_rates[2136] - rev_rates[995]);
  //sp 117
  sp_rates[116] += (fwd_rates[2136] - rev_rates[995]);
  //sp 1
  (*dy_N) += (fwd_rates[2136] - rev_rates[995]);

  //rxn 2137
  //sp 10
  sp_rates[9] -= (fwd_rates[2137] - rev_rates[996]);
  //sp 116
  sp_rates[115] += (fwd_rates[2137] - rev_rates[996]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2137] - rev_rates[996]);
  //sp 8
  sp_rates[7] += (fwd_rates[2137] - rev_rates[996]);

  //rxn 2138
  //sp 114
  sp_rates[113] -= (fwd_rates[2138] - rev_rates[997]) * pres_mod[53];
  //sp 116
  sp_rates[115] += (fwd_rates[2138] - rev_rates[997]) * pres_mod[53];
  //sp 6
  sp_rates[5] -= (fwd_rates[2138] - rev_rates[997]) * pres_mod[53];

  //rxn 2139
  //sp 114
  sp_rates[113] -= (fwd_rates[2139] - rev_rates[998]);
  //sp 119
  sp_rates[118] += (fwd_rates[2139] - rev_rates[998]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2139] - rev_rates[998]);

  //rxn 2140
  //sp 26
  sp_rates[25] -= (fwd_rates[2140] - rev_rates[999]);
  //sp 11
  sp_rates[10] += (fwd_rates[2140] - rev_rates[999]);
  //sp 117
  sp_rates[116] += (fwd_rates[2140] - rev_rates[999]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2140] - rev_rates[999]);

  //rxn 2141
  //sp 114
  sp_rates[113] -= (fwd_rates[2141] - rev_rates[1000]) * pres_mod[54];
  //sp 4
  sp_rates[3] -= (fwd_rates[2141] - rev_rates[1000]) * pres_mod[54];
  //sp 117
  sp_rates[116] += (fwd_rates[2141] - rev_rates[1000]) * pres_mod[54];

  //rxn 2142
  //sp 114
  sp_rates[113] += (fwd_rates[2142] - rev_rates[1001]);
  //sp 3
  sp_rates[2] += (fwd_rates[2142] - rev_rates[1001]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2142] - rev_rates[1001]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2142] - rev_rates[1001]);

  //rxn 2143
  //sp 114
  sp_rates[113] += (fwd_rates[2143] - rev_rates[1002]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2143] - rev_rates[1002]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2143] - rev_rates[1002]);
  //sp 8
  sp_rates[7] += (fwd_rates[2143] - rev_rates[1002]);

  //rxn 2144
  //sp 4
  sp_rates[3] += (fwd_rates[2144] - rev_rates[1003]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2144] - rev_rates[1003]);
  //sp 119
  sp_rates[118] += (fwd_rates[2144] - rev_rates[1003]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2144] - rev_rates[1003]);

  //rxn 2145
  //sp 114
  sp_rates[113] += (fwd_rates[2145] - rev_rates[1004]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2145] - rev_rates[1004]);
  //sp 7
  sp_rates[6] += (fwd_rates[2145] - rev_rates[1004]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2145] - rev_rates[1004]);

  //rxn 2146
  //sp 117
  sp_rates[116] -= (fwd_rates[2146] - rev_rates[1005]);
  //sp 10
  sp_rates[9] += (fwd_rates[2146] - rev_rates[1005]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2146] - rev_rates[1005]);
  //sp 114
  sp_rates[113] += (fwd_rates[2146] - rev_rates[1005]);

  //rxn 2147
  //sp 133
  sp_rates[132] -= (fwd_rates[2147] - rev_rates[1006]);
  //sp 127
  sp_rates[126] += (fwd_rates[2147] - rev_rates[1006]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2147] - rev_rates[1006]);
  //sp 114
  sp_rates[113] += (fwd_rates[2147] - rev_rates[1006]);

  //rxn 2148
  //sp 124
  sp_rates[123] += (fwd_rates[2148] - rev_rates[1007]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2148] - rev_rates[1007]);
  //sp 4
  sp_rates[3] += (fwd_rates[2148] - rev_rates[1007]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2148] - rev_rates[1007]);

  //rxn 2149
  //sp 114
  sp_rates[113] -= (fwd_rates[2149] - rev_rates[1008]);
  //sp 115
  sp_rates[114] += (fwd_rates[2149] - rev_rates[1008]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2149] - rev_rates[1008]);
  //sp 8
  sp_rates[7] += (fwd_rates[2149] - rev_rates[1008]);

  //rxn 2150
  //sp 114
  sp_rates[113] += (fwd_rates[2150] - rev_rates[1009]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2150] - rev_rates[1009]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2150] - rev_rates[1009]);
  //sp 119
  sp_rates[118] += (fwd_rates[2150] - rev_rates[1009]);

  //rxn 2151
  //sp 114
  sp_rates[113] += (fwd_rates[2151] - rev_rates[1010]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2151] - rev_rates[1010]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2151] - rev_rates[1010]);
  //sp 118
  sp_rates[117] += (fwd_rates[2151] - rev_rates[1010]);

  //rxn 2152
  //sp 115
  sp_rates[114] += (fwd_rates[2152] - rev_rates[1011]);
  //sp 117
  sp_rates[116] -= 2.0 * (fwd_rates[2152] - rev_rates[1011]);
  //sp 7
  sp_rates[6] += (fwd_rates[2152] - rev_rates[1011]);

  //rxn 2153
  //sp 114
  sp_rates[113] += (fwd_rates[2153] - rev_rates[1012]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2153] - rev_rates[1012]);
  //sp 14
  sp_rates[13] += (fwd_rates[2153] - rev_rates[1012]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2153] - rev_rates[1012]);

  //rxn 2154
  //sp 116
  sp_rates[115] += (fwd_rates[2154] - rev_rates[1013]);
  //sp 3
  sp_rates[2] += (fwd_rates[2154] - rev_rates[1013]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2154] - rev_rates[1013]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2154] - rev_rates[1013]);

  //rxn 2155
  //sp 116
  sp_rates[115] += (fwd_rates[2155] - rev_rates[1014]);
  //sp 3
  sp_rates[2] += (fwd_rates[2155] - rev_rates[1014]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2155] - rev_rates[1014]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2155] - rev_rates[1014]);

  //rxn 2156
  //sp 114
  sp_rates[113] += (fwd_rates[2156] - rev_rates[1015]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2156] - rev_rates[1015]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2156] - rev_rates[1015]);
  //sp 7
  sp_rates[6] += (fwd_rates[2156] - rev_rates[1015]);

  //rxn 2157
  //sp 119
  sp_rates[118] -= (fwd_rates[2157] - rev_rates[1016]);
  //sp 116
  sp_rates[115] += (fwd_rates[2157] - rev_rates[1016]);
  //sp 14
  sp_rates[13] += (fwd_rates[2157] - rev_rates[1016]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2157] - rev_rates[1016]);

  //rxn 2158
  //sp 119
  sp_rates[118] -= (fwd_rates[2158] - rev_rates[1017]);
  //sp 116
  sp_rates[115] += (fwd_rates[2158] - rev_rates[1017]);
  //sp 14
  sp_rates[13] += (fwd_rates[2158] - rev_rates[1017]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2158] - rev_rates[1017]);

  //rxn 2159
  //sp 14
  sp_rates[13] += (fwd_rates[2159] - rev_rates[1018]);
  //sp 116
  sp_rates[115] += (fwd_rates[2159] - rev_rates[1018]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2159] - rev_rates[1018]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2159] - rev_rates[1018]);

  //rxn 2160
  //sp 14
  sp_rates[13] += (fwd_rates[2160] - rev_rates[1019]);
  //sp 116
  sp_rates[115] += (fwd_rates[2160] - rev_rates[1019]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2160] - rev_rates[1019]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2160] - rev_rates[1019]);

  //rxn 2161
  //sp 114
  sp_rates[113] += (fwd_rates[2161] - rev_rates[1020]);
  //sp 22
  sp_rates[21] += (fwd_rates[2161] - rev_rates[1020]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2161] - rev_rates[1020]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2161] - rev_rates[1020]);

  //rxn 2162
  //sp 141
  sp_rates[140] = (fwd_rates[2162] - rev_rates[1021]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2162] - rev_rates[1021]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2162] - rev_rates[1021]);
  //sp 7
  sp_rates[6] += (fwd_rates[2162] - rev_rates[1021]);

  //rxn 2163
  //sp 114
  sp_rates[113] += (fwd_rates[2163] - rev_rates[1022]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2163] - rev_rates[1022]);
  //sp 22
  sp_rates[21] += (fwd_rates[2163] - rev_rates[1022]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2163] - rev_rates[1022]);

  //rxn 2164
  //sp 119
  sp_rates[118] -= (fwd_rates[2164] - rev_rates[1023]);
  //sp 141
  sp_rates[140] += (fwd_rates[2164] - rev_rates[1023]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2164] - rev_rates[1023]);
  //sp 7
  sp_rates[6] += (fwd_rates[2164] - rev_rates[1023]);

  //rxn 2165
  //sp 114
  sp_rates[113] += (fwd_rates[2165] - rev_rates[1024]);
  //sp 22
  sp_rates[21] += (fwd_rates[2165] - rev_rates[1024]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2165] - rev_rates[1024]);
  //sp 14
  sp_rates[13] -= (fwd_rates[2165] - rev_rates[1024]);

  //rxn 2166
  //sp 116
  sp_rates[115] -= (fwd_rates[2166] - rev_rates[1025]);
  //sp 141
  sp_rates[140] += (fwd_rates[2166] - rev_rates[1025]);
  //sp 14
  sp_rates[13] -= (fwd_rates[2166] - rev_rates[1025]);
  //sp 7
  sp_rates[6] += (fwd_rates[2166] - rev_rates[1025]);

  //rxn 2167
  //sp 119
  sp_rates[118] += (fwd_rates[2167] - rev_rates[1026]);
  //sp 142
  sp_rates[141] = -(fwd_rates[2167] - rev_rates[1026]);
  //sp 15
  sp_rates[14] += (fwd_rates[2167] - rev_rates[1026]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2167] - rev_rates[1026]);

  //rxn 2168
  //sp 114
  sp_rates[113] += (fwd_rates[2168] - rev_rates[1027]);
  //sp 22
  sp_rates[21] += (fwd_rates[2168] - rev_rates[1027]);
  //sp 142
  sp_rates[141] -= (fwd_rates[2168] - rev_rates[1027]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2168] - rev_rates[1027]);

  //rxn 2169
  //sp 114
  sp_rates[113] += (fwd_rates[2169] - rev_rates[1028]);
  //sp 141
  sp_rates[140] -= (fwd_rates[2169] - rev_rates[1028]);
  //sp 22
  sp_rates[21] += (fwd_rates[2169] - rev_rates[1028]);
  //sp 7
  sp_rates[6] -= (fwd_rates[2169] - rev_rates[1028]);

  //rxn 2170
  //sp 127
  sp_rates[126] -= (fwd_rates[2170] - rev_rates[1029]);
  //sp 133
  sp_rates[132] += (fwd_rates[2170] - rev_rates[1029]);
  //sp 14
  sp_rates[13] += (fwd_rates[2170] - rev_rates[1029]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2170] - rev_rates[1029]);

  //rxn 2171
  //sp 116
  sp_rates[115] += (fwd_rates[2171] - rev_rates[1030]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2171] - rev_rates[1030]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2171] - rev_rates[1030]);
  //sp 8
  sp_rates[7] += (fwd_rates[2171] - rev_rates[1030]);

  //rxn 2172
  //sp 116
  sp_rates[115] += (fwd_rates[2172] - rev_rates[1031]);
  //sp 7
  sp_rates[6] += (fwd_rates[2172] - rev_rates[1031]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2172] - rev_rates[1031]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2172] - rev_rates[1031]);

  //rxn 2173
  //sp 131
  sp_rates[130] -= (fwd_rates[2173] - rev_rates[1032]);
  //sp 116
  sp_rates[115] += (fwd_rates[2173] - rev_rates[1032]);
  //sp 133
  sp_rates[132] += (fwd_rates[2173] - rev_rates[1032]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2173] - rev_rates[1032]);

  //rxn 2174
  //sp 124
  sp_rates[123] += (fwd_rates[2174] - rev_rates[1033]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2174] - rev_rates[1033]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2174] - rev_rates[1033]);
  //sp 8
  sp_rates[7] += (fwd_rates[2174] - rev_rates[1033]);

  //rxn 2175
  //sp 127
  sp_rates[126] += (fwd_rates[2175] - rev_rates[1034]);
  //sp 116
  sp_rates[115] += (fwd_rates[2175] - rev_rates[1034]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2175] - rev_rates[1034]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2175] - rev_rates[1034]);

  //rxn 2176
  //sp 114
  sp_rates[113] += (fwd_rates[2176] - rev_rates[1035]);
  //sp 116
  sp_rates[115] += (fwd_rates[2176] - rev_rates[1035]);
  //sp 119
  sp_rates[118] -= 2.0 * (fwd_rates[2176] - rev_rates[1035]);
  //sp 7
  sp_rates[6] += (fwd_rates[2176] - rev_rates[1035]);

  //rxn 2177
  //sp 126
  sp_rates[125] += (fwd_rates[2177] - rev_rates[1036]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2177] - rev_rates[1036]);

  //rxn 2178
  //sp 4
  sp_rates[3] += (fwd_rates[2178] - rev_rates[1037]) * pres_mod[55];
  //sp 117
  sp_rates[116] += (fwd_rates[2178] - rev_rates[1037]) * pres_mod[55];
  //sp 134
  sp_rates[133] -= (fwd_rates[2178] - rev_rates[1037]) * pres_mod[55];

  //rxn 2179
  //sp 3
  sp_rates[2] += (fwd_rates[2179] - rev_rates[1038]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2179] - rev_rates[1038]);
  //sp 117
  sp_rates[116] += (fwd_rates[2179] - rev_rates[1038]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2179] - rev_rates[1038]);

  //rxn 2180
  //sp 4
  sp_rates[3] -= (fwd_rates[2180] - rev_rates[1039]);
  //sp 133
  sp_rates[132] += (fwd_rates[2180] - rev_rates[1039]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2180] - rev_rates[1039]);
  //sp 8
  sp_rates[7] += (fwd_rates[2180] - rev_rates[1039]);

  //rxn 2181
  //sp 117
  sp_rates[116] += (fwd_rates[2181] - rev_rates[1040]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2181] - rev_rates[1040]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2181] - rev_rates[1040]);
  //sp 8
  sp_rates[7] += (fwd_rates[2181] - rev_rates[1040]);

  //rxn 2182
  //sp 117
  sp_rates[116] += (fwd_rates[2182] - rev_rates[1041]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2182] - rev_rates[1041]);
  //sp 7
  sp_rates[6] += (fwd_rates[2182] - rev_rates[1041]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2182] - rev_rates[1041]);

  //rxn 2183
  //sp 133
  sp_rates[132] += (fwd_rates[2183] - rev_rates[1042]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2183] - rev_rates[1042]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2183] - rev_rates[1042]);
  //sp 23
  sp_rates[22] += (fwd_rates[2183] - rev_rates[1042]);

  //rxn 2184
  //sp 14
  sp_rates[13] += (fwd_rates[2184] - rev_rates[1043]);
  //sp 117
  sp_rates[116] += (fwd_rates[2184] - rev_rates[1043]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2184] - rev_rates[1043]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2184] - rev_rates[1043]);

  //rxn 2185
  //sp 114
  sp_rates[113] -= (fwd_rates[2185] - rev_rates[1044]);
  //sp 117
  sp_rates[116] += 2.0 * (fwd_rates[2185] - rev_rates[1044]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2185] - rev_rates[1044]);

  //rxn 2186
  //sp 116
  sp_rates[115] -= (fwd_rates[2186] - rev_rates[1045]);
  //sp 117
  sp_rates[116] += (fwd_rates[2186] - rev_rates[1045]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2186] - rev_rates[1045]);
  //sp 119
  sp_rates[118] += (fwd_rates[2186] - rev_rates[1045]);

  //rxn 2187
  //sp 117
  sp_rates[116] += (fwd_rates[2187] - rev_rates[1046]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2187] - rev_rates[1046]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2187] - rev_rates[1046]);
  //sp 127
  sp_rates[126] += (fwd_rates[2187] - rev_rates[1046]);

  //rxn 2188
  //sp 117
  sp_rates[116] += (fwd_rates[2188] - rev_rates[1047]);
  //sp 10
  sp_rates[9] += (fwd_rates[2188] - rev_rates[1047]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2188] - rev_rates[1047]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2188] - rev_rates[1047]);

  //rxn 2189
  //sp 9
  sp_rates[8] += (fwd_rates[2189] - rev_rates[1048]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2189] - rev_rates[1048]);
  //sp 117
  sp_rates[116] += (fwd_rates[2189] - rev_rates[1048]);
  //sp 134
  sp_rates[133] -= (fwd_rates[2189] - rev_rates[1048]);

  //rxn 2190
  //sp 116
  sp_rates[115] += (fwd_rates[2190] - rev_rates[1049]);
  //sp 130
  sp_rates[129] = -(fwd_rates[2190] - rev_rates[1049]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2190] - rev_rates[1049]);
  //sp 8
  sp_rates[7] += (fwd_rates[2190] - rev_rates[1049]);

  //rxn 2191
  //sp 130
  sp_rates[129] -= (fwd_rates[2191] - rev_rates[1050]);
  //sp 116
  sp_rates[115] += (fwd_rates[2191] - rev_rates[1050]);
  //sp 5
  sp_rates[4] += (fwd_rates[2191] - rev_rates[1050]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2191] - rev_rates[1050]);

  //rxn 2192
  //sp 130
  sp_rates[129] -= (fwd_rates[2192] - rev_rates[1051]);
  //sp 116
  sp_rates[115] += (fwd_rates[2192] - rev_rates[1051]);
  //sp 10
  sp_rates[9] += (fwd_rates[2192] - rev_rates[1051]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2192] - rev_rates[1051]);

  //rxn 2193
  //sp 130
  sp_rates[129] -= (fwd_rates[2193] - rev_rates[1052]);
  //sp 5
  sp_rates[4] += (fwd_rates[2193] - rev_rates[1052]);
  //sp 8
  sp_rates[7] += (fwd_rates[2193] - rev_rates[1052]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2193] - rev_rates[1052]);
  //sp 116
  sp_rates[115] += (fwd_rates[2193] - rev_rates[1052]);

  //rxn 2194
  //sp 130
  sp_rates[129] -= (fwd_rates[2194] - rev_rates[1053]);
  //sp 5
  sp_rates[4] += (fwd_rates[2194] - rev_rates[1053]);
  //sp 114
  sp_rates[113] += (fwd_rates[2194] - rev_rates[1053]);

  //rxn 2195
  //sp 130
  sp_rates[129] -= (fwd_rates[2195] - rev_rates[1054]);
  //sp 5
  sp_rates[4] += (fwd_rates[2195] - rev_rates[1054]);
  //sp 114
  sp_rates[113] += (fwd_rates[2195] - rev_rates[1054]);

  //rxn 2196
  //sp 130
  sp_rates[129] -= 2.0 * (fwd_rates[2196] - rev_rates[1055]);
  //sp 116
  sp_rates[115] += 2.0 * (fwd_rates[2196] - rev_rates[1055]);
  //sp 5
  sp_rates[4] += (fwd_rates[2196] - rev_rates[1055]);

  //rxn 2197
  //sp 116
  sp_rates[115] -= (fwd_rates[2197] - rev_rates[1056]);
  //sp 114
  sp_rates[113] += (fwd_rates[2197] - rev_rates[1056]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2197] - rev_rates[1056]);
  //sp 8
  sp_rates[7] += (fwd_rates[2197] - rev_rates[1056]);

  //rxn 2198
  //sp 114
  sp_rates[113] += (fwd_rates[2198] - rev_rates[1057]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2198] - rev_rates[1057]);
  //sp 5
  sp_rates[4] += (fwd_rates[2198] - rev_rates[1057]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2198] - rev_rates[1057]);

  //rxn 2199
  //sp 130
  sp_rates[129] += (fwd_rates[2199] - rev_rates[1058]) * pres_mod[56];
  //sp 116
  sp_rates[115] -= (fwd_rates[2199] - rev_rates[1058]) * pres_mod[56];
  //sp 6
  sp_rates[5] -= (fwd_rates[2199] - rev_rates[1058]) * pres_mod[56];

  //rxn 2200
  //sp 10
  sp_rates[9] -= (fwd_rates[2200] - rev_rates[1059]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2200] - rev_rates[1059]);
  //sp 5
  sp_rates[4] += (fwd_rates[2200] - rev_rates[1059]);
  //sp 119
  sp_rates[118] += (fwd_rates[2200] - rev_rates[1059]);

  //rxn 2201
  //sp 10
  sp_rates[9] -= (fwd_rates[2201] - rev_rates[1060]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2201] - rev_rates[1060]);
  //sp 5
  sp_rates[4] += (fwd_rates[2201] - rev_rates[1060]);
  //sp 118
  sp_rates[117] += (fwd_rates[2201] - rev_rates[1060]);

  //rxn 2202
  //sp 114
  sp_rates[113] += 2.0 * (fwd_rates[2202] - rev_rates[1061]);
  //sp 116
  sp_rates[115] -= 2.0 * (fwd_rates[2202] - rev_rates[1061]);
  //sp 5
  sp_rates[4] += (fwd_rates[2202] - rev_rates[1061]);

  //rxn 2203
  //sp 114
  sp_rates[113] += (fwd_rates[2203] - rev_rates[1062]);
  //sp 130
  sp_rates[129] += (fwd_rates[2203] - rev_rates[1062]);
  //sp 116
  sp_rates[115] -= 2.0 * (fwd_rates[2203] - rev_rates[1062]);

  //rxn 2204
  //sp 12
  sp_rates[11] += (fwd_rates[2204] - rev_rates[1063]);
  //sp 114
  sp_rates[113] += (fwd_rates[2204] - rev_rates[1063]);
  //sp 11
  sp_rates[10] -= (fwd_rates[2204] - rev_rates[1063]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2204] - rev_rates[1063]);

  //rxn 2205
  //sp 4
  sp_rates[3] += (fwd_rates[2205] - rev_rates[1064]);
  //sp 12
  sp_rates[11] += (fwd_rates[2205] - rev_rates[1064]);
  //sp 114
  sp_rates[113] += (fwd_rates[2205] - rev_rates[1064]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2205] - rev_rates[1064]);
  //sp 26
  sp_rates[25] -= (fwd_rates[2205] - rev_rates[1064]);

  //rxn 2206
  //sp 26
  sp_rates[25] -= (fwd_rates[2206] - rev_rates[1065]);
  //sp 11
  sp_rates[10] += (fwd_rates[2206] - rev_rates[1065]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2206] - rev_rates[1065]);
  //sp 119
  sp_rates[118] += (fwd_rates[2206] - rev_rates[1065]);

  //rxn 2207
  //sp 8
  sp_rates[7] += (fwd_rates[2207] - rev_rates[1066]);
  //sp 11
  sp_rates[10] += (fwd_rates[2207] - rev_rates[1066]);
  //sp 114
  sp_rates[113] += (fwd_rates[2207] - rev_rates[1066]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2207] - rev_rates[1066]);
  //sp 26
  sp_rates[25] -= (fwd_rates[2207] - rev_rates[1066]);

  //rxn 2208
  //sp 4
  sp_rates[3] += (fwd_rates[2208] - rev_rates[1067]);
  //sp 149
  sp_rates[148] = (fwd_rates[2208] - rev_rates[1067]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2208] - rev_rates[1067]);
  //sp 136
  sp_rates[135] = -(fwd_rates[2208] - rev_rates[1067]);

  //rxn 2209
  //sp 11
  sp_rates[10] += (fwd_rates[2209] - rev_rates[1068]);
  //sp 131
  sp_rates[130] += (fwd_rates[2209] - rev_rates[1068]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2209] - rev_rates[1068]);
  //sp 136
  sp_rates[135] -= (fwd_rates[2209] - rev_rates[1068]);

  //rxn 2210
  //sp 147
  sp_rates[146] = (fwd_rates[2210] - rev_rates[1069]);
  //sp 8
  sp_rates[7] += (fwd_rates[2210] - rev_rates[1069]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2210] - rev_rates[1069]);
  //sp 136
  sp_rates[135] -= (fwd_rates[2210] - rev_rates[1069]);

  //rxn 2211
  //sp 136
  sp_rates[135] -= (fwd_rates[2211] - rev_rates[1070]);
  //sp 147
  sp_rates[146] += (fwd_rates[2211] - rev_rates[1070]);
  //sp 7
  sp_rates[6] += (fwd_rates[2211] - rev_rates[1070]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2211] - rev_rates[1070]);

  //rxn 2212
  //sp 136
  sp_rates[135] -= (fwd_rates[2212] - rev_rates[1071]);
  //sp 138
  sp_rates[137] = (fwd_rates[2212] - rev_rates[1071]);
  //sp 4
  sp_rates[3] += (fwd_rates[2212] - rev_rates[1071]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2212] - rev_rates[1071]);

  //rxn 2213
  //sp 136
  sp_rates[135] -= (fwd_rates[2213] - rev_rates[1072]);
  //sp 11
  sp_rates[10] += (fwd_rates[2213] - rev_rates[1072]);
  //sp 133
  sp_rates[132] += (fwd_rates[2213] - rev_rates[1072]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2213] - rev_rates[1072]);

  //rxn 2214
  //sp 10
  sp_rates[9] += (fwd_rates[2214] - rev_rates[1073]);
  //sp 147
  sp_rates[146] += (fwd_rates[2214] - rev_rates[1073]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2214] - rev_rates[1073]);
  //sp 136
  sp_rates[135] -= (fwd_rates[2214] - rev_rates[1073]);

  //rxn 2215
  //sp 136
  sp_rates[135] += (fwd_rates[2215] - rev_rates[1074]);
  //sp 3
  sp_rates[2] -= (fwd_rates[2215] - rev_rates[1074]);
  //sp 4
  sp_rates[3] += (fwd_rates[2215] - rev_rates[1074]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2215] - rev_rates[1074]);

  //rxn 2216
  //sp 129
  sp_rates[128] += (fwd_rates[2216] - rev_rates[1075]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2216] - rev_rates[1075]);
  //sp 11
  sp_rates[10] += (fwd_rates[2216] - rev_rates[1075]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2216] - rev_rates[1075]);

  //rxn 2217
  //sp 147
  sp_rates[146] -= (fwd_rates[2217] - rev_rates[1076]);
  //sp 4
  sp_rates[3] += (fwd_rates[2217] - rev_rates[1076]);
  //sp 149
  sp_rates[148] += (fwd_rates[2217] - rev_rates[1076]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2217] - rev_rates[1076]);

  //rxn 2218
  //sp 114
  sp_rates[113] += (fwd_rates[2218] - rev_rates[1077]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2218] - rev_rates[1077]);
  //sp 11
  sp_rates[10] += (fwd_rates[2218] - rev_rates[1077]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2218] - rev_rates[1077]);

  //rxn 2219
  //sp 149
  sp_rates[148] += (fwd_rates[2219] - rev_rates[1078]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2219] - rev_rates[1078]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2219] - rev_rates[1078]);
  //sp 6
  sp_rates[5] += (fwd_rates[2219] - rev_rates[1078]);

  //rxn 2220
  //sp 129
  sp_rates[128] += (fwd_rates[2220] - rev_rates[1079]) * pres_mod[57];
  //sp 11
  sp_rates[10] += (fwd_rates[2220] - rev_rates[1079]) * pres_mod[57];
  //sp 149
  sp_rates[148] -= (fwd_rates[2220] - rev_rates[1079]) * pres_mod[57];

  //rxn 2221
  //sp 11
  sp_rates[10] += (fwd_rates[2221] - rev_rates[1080]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2221] - rev_rates[1080]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2221] - rev_rates[1080]);
  //sp 131
  sp_rates[130] += (fwd_rates[2221] - rev_rates[1080]);

  //rxn 2222
  //sp 114
  sp_rates[113] += (fwd_rates[2222] - rev_rates[1081]);
  //sp 11
  sp_rates[10] += (fwd_rates[2222] - rev_rates[1081]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2222] - rev_rates[1081]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2222] - rev_rates[1081]);

  //rxn 2223
  //sp 129
  sp_rates[128] -= (fwd_rates[2223] - rev_rates[1082]);
  //sp 11
  sp_rates[10] += (fwd_rates[2223] - rev_rates[1082]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2223] - rev_rates[1082]);
  //sp 1
  (*dy_N) += (fwd_rates[2223] - rev_rates[1082]);

  //rxn 2224
  //sp 149
  sp_rates[148] -= (fwd_rates[2224] - rev_rates[1083]);
  //sp 12
  sp_rates[11] += (fwd_rates[2224] - rev_rates[1083]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2224] - rev_rates[1083]);
  //sp 114
  sp_rates[113] += (fwd_rates[2224] - rev_rates[1083]);

  //rxn 2225
  //sp 138
  sp_rates[137] += (fwd_rates[2225] - rev_rates[1084]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2225] - rev_rates[1084]);
  //sp 6
  sp_rates[5] += (fwd_rates[2225] - rev_rates[1084]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2225] - rev_rates[1084]);

  //rxn 2226
  //sp 4
  sp_rates[3] += (fwd_rates[2226] - rev_rates[1085]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2226] - rev_rates[1085]);
  //sp 11
  sp_rates[10] += (fwd_rates[2226] - rev_rates[1085]);
  //sp 114
  sp_rates[113] += (fwd_rates[2226] - rev_rates[1085]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2226] - rev_rates[1085]);

  //rxn 2227
  //sp 26
  sp_rates[25] -= (fwd_rates[2227] - rev_rates[1086]);
  //sp 11
  sp_rates[10] += (fwd_rates[2227] - rev_rates[1086]);
  //sp 138
  sp_rates[137] += (fwd_rates[2227] - rev_rates[1086]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2227] - rev_rates[1086]);

  //rxn 2228
  //sp 25
  sp_rates[24] -= (fwd_rates[2228] - rev_rates[1087]);
  //sp 26
  sp_rates[25] += (fwd_rates[2228] - rev_rates[1087]);
  //sp 138
  sp_rates[137] += (fwd_rates[2228] - rev_rates[1087]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2228] - rev_rates[1087]);

  //rxn 2229
  //sp 114
  sp_rates[113] += 2.0 * (fwd_rates[2229] - rev_rates[1088]);
  //sp 11
  sp_rates[10] += (fwd_rates[2229] - rev_rates[1088]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2229] - rev_rates[1088]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2229] - rev_rates[1088]);

  //rxn 2230
  //sp 12
  sp_rates[11] += (fwd_rates[2230] - rev_rates[1089]);
  //sp 115
  sp_rates[114] += (fwd_rates[2230] - rev_rates[1089]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2230] - rev_rates[1089]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2230] - rev_rates[1089]);

  //rxn 2231
  //sp 149
  sp_rates[148] -= (fwd_rates[2231] - rev_rates[1090]);
  //sp 138
  sp_rates[137] += (fwd_rates[2231] - rev_rates[1090]);
  //sp 117
  sp_rates[116] -= (fwd_rates[2231] - rev_rates[1090]);
  //sp 114
  sp_rates[113] += (fwd_rates[2231] - rev_rates[1090]);

  //rxn 2232
  //sp 138
  sp_rates[137] += (fwd_rates[2232] - rev_rates[1091]);
  //sp 116
  sp_rates[115] += (fwd_rates[2232] - rev_rates[1091]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2232] - rev_rates[1091]);
  //sp 119
  sp_rates[118] -= (fwd_rates[2232] - rev_rates[1091]);

  //rxn 2233
  //sp 11
  sp_rates[10] += 2.0 * (fwd_rates[2233] - rev_rates[1092]);
  //sp 149
  sp_rates[148] -= 2.0 * (fwd_rates[2233] - rev_rates[1092]);
  //sp 1
  (*dy_N) += (fwd_rates[2233] - rev_rates[1092]);

  //rxn 2234
  //sp 133
  sp_rates[132] += (fwd_rates[2234] - rev_rates[1093]);
  //sp 138
  sp_rates[137] += (fwd_rates[2234] - rev_rates[1093]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2234] - rev_rates[1093]);
  //sp 127
  sp_rates[126] -= (fwd_rates[2234] - rev_rates[1093]);

  //rxn 2235
  //sp 138
  sp_rates[137] += (fwd_rates[2235] - rev_rates[1094]);
  //sp 140
  sp_rates[139] = -(fwd_rates[2235] - rev_rates[1094]);

  //rxn 2236
  //sp 136
  sp_rates[135] += (fwd_rates[2236] - rev_rates[1095]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2236] - rev_rates[1095]);
  //sp 140
  sp_rates[139] -= (fwd_rates[2236] - rev_rates[1095]);
  //sp 8
  sp_rates[7] += (fwd_rates[2236] - rev_rates[1095]);

  //rxn 2237
  //sp 133
  sp_rates[132] += (fwd_rates[2237] - rev_rates[1096]);
  //sp 11
  sp_rates[10] += (fwd_rates[2237] - rev_rates[1096]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2237] - rev_rates[1096]);
  //sp 140
  sp_rates[139] -= (fwd_rates[2237] - rev_rates[1096]);

  //rxn 2238
  //sp 149
  sp_rates[148] += (fwd_rates[2238] - rev_rates[1097]);
  //sp 3
  sp_rates[2] += (fwd_rates[2238] - rev_rates[1097]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2238] - rev_rates[1097]);
  //sp 140
  sp_rates[139] -= (fwd_rates[2238] - rev_rates[1097]);

  //rxn 2239
  //sp 140
  sp_rates[139] -= (fwd_rates[2239] - rev_rates[1098]);
  //sp 149
  sp_rates[148] += (fwd_rates[2239] - rev_rates[1098]);
  //sp 7
  sp_rates[6] += (fwd_rates[2239] - rev_rates[1098]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2239] - rev_rates[1098]);

  //rxn 2240
  //sp 140
  sp_rates[139] -= (fwd_rates[2240] - rev_rates[1099]);
  //sp 149
  sp_rates[148] += (fwd_rates[2240] - rev_rates[1099]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2240] - rev_rates[1099]);
  //sp 8
  sp_rates[7] += (fwd_rates[2240] - rev_rates[1099]);

  //rxn 2241
  //sp 136
  sp_rates[135] += (fwd_rates[2241] - rev_rates[1100]);
  //sp 139
  sp_rates[138] = -(fwd_rates[2241] - rev_rates[1100]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2241] - rev_rates[1100]);
  //sp 8
  sp_rates[7] += (fwd_rates[2241] - rev_rates[1100]);

  //rxn 2242
  //sp 138
  sp_rates[137] += (fwd_rates[2242] - rev_rates[1101]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2242] - rev_rates[1101]);

  //rxn 2243
  //sp 11
  sp_rates[10] += (fwd_rates[2243] - rev_rates[1102]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2243] - rev_rates[1102]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2243] - rev_rates[1102]);
  //sp 133
  sp_rates[132] += (fwd_rates[2243] - rev_rates[1102]);

  //rxn 2244
  //sp 139
  sp_rates[138] -= (fwd_rates[2244] - rev_rates[1103]);
  //sp 140
  sp_rates[139] += (fwd_rates[2244] - rev_rates[1103]);

  //rxn 2245
  //sp 139
  sp_rates[138] -= (fwd_rates[2245] - rev_rates[1104]);
  //sp 149
  sp_rates[148] += (fwd_rates[2245] - rev_rates[1104]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2245] - rev_rates[1104]);
  //sp 8
  sp_rates[7] += (fwd_rates[2245] - rev_rates[1104]);

  //rxn 2246
  //sp 26
  sp_rates[25] += (fwd_rates[2246] - rev_rates[1105]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2246] - rev_rates[1105]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2246] - rev_rates[1105]);
  //sp 114
  sp_rates[113] += (fwd_rates[2246] - rev_rates[1105]);

  //rxn 2247
  //sp 25
  sp_rates[24] += (fwd_rates[2247] - rev_rates[1106]);
  //sp 114
  sp_rates[113] += (fwd_rates[2247] - rev_rates[1106]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2247] - rev_rates[1106]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2247] - rev_rates[1106]);

  //rxn 2248
  //sp 3
  sp_rates[2] += (fwd_rates[2248] - rev_rates[1107]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2248] - rev_rates[1107]);
  //sp 11
  sp_rates[10] += (fwd_rates[2248] - rev_rates[1107]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2248] - rev_rates[1107]);
  //sp 114
  sp_rates[113] += (fwd_rates[2248] - rev_rates[1107]);

  //rxn 2249
  //sp 139
  sp_rates[138] -= fwd_rates[2249];
  //sp 4
  sp_rates[3] += fwd_rates[2249];
  //sp 149
  sp_rates[148] += fwd_rates[2249];

  //rxn 2250
  //sp 139
  sp_rates[138] -= (fwd_rates[2250] - rev_rates[1108]);
  //sp 149
  sp_rates[148] += (fwd_rates[2250] - rev_rates[1108]);
  //sp 7
  sp_rates[6] += (fwd_rates[2250] - rev_rates[1108]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2250] - rev_rates[1108]);

  //rxn 2251
  //sp 26
  sp_rates[25] += (fwd_rates[2251] - rev_rates[1109]);
  //sp 139
  sp_rates[138] -= (fwd_rates[2251] - rev_rates[1109]);
  //sp 117
  sp_rates[116] += (fwd_rates[2251] - rev_rates[1109]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2251] - rev_rates[1109]);

  //rxn 2252
  //sp 138
  sp_rates[137] -= (fwd_rates[2252] - rev_rates[1110]) * pres_mod[58];
  //sp 11
  sp_rates[10] += (fwd_rates[2252] - rev_rates[1110]) * pres_mod[58];
  //sp 131
  sp_rates[130] += (fwd_rates[2252] - rev_rates[1110]) * pres_mod[58];

  //rxn 2253
  //sp 138
  sp_rates[137] -= (fwd_rates[2253] - rev_rates[1111]);
  //sp 11
  sp_rates[10] += (fwd_rates[2253] - rev_rates[1111]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2253] - rev_rates[1111]);
  //sp 133
  sp_rates[132] += (fwd_rates[2253] - rev_rates[1111]);

  //rxn 2254
  //sp 138
  sp_rates[137] -= (fwd_rates[2254] - rev_rates[1112]);
  //sp 3
  sp_rates[2] += (fwd_rates[2254] - rev_rates[1112]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2254] - rev_rates[1112]);
  //sp 149
  sp_rates[148] += (fwd_rates[2254] - rev_rates[1112]);

  //rxn 2255
  //sp 138
  sp_rates[137] -= (fwd_rates[2255] - rev_rates[1113]);
  //sp 131
  sp_rates[130] += (fwd_rates[2255] - rev_rates[1113]);
  //sp 12
  sp_rates[11] += (fwd_rates[2255] - rev_rates[1113]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2255] - rev_rates[1113]);

  //rxn 2256
  //sp 138
  sp_rates[137] -= (fwd_rates[2256] - rev_rates[1114]);
  //sp 11
  sp_rates[10] += (fwd_rates[2256] - rev_rates[1114]);
  //sp 117
  sp_rates[116] += (fwd_rates[2256] - rev_rates[1114]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2256] - rev_rates[1114]);

  //rxn 2257
  //sp 138
  sp_rates[137] -= (fwd_rates[2257] - rev_rates[1115]);
  //sp 149
  sp_rates[148] += (fwd_rates[2257] - rev_rates[1115]);
  //sp 7
  sp_rates[6] += (fwd_rates[2257] - rev_rates[1115]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2257] - rev_rates[1115]);

  //rxn 2258
  //sp 9
  sp_rates[8] += (fwd_rates[2258] - rev_rates[1116]);
  //sp 10
  sp_rates[9] -= (fwd_rates[2258] - rev_rates[1116]);
  //sp 138
  sp_rates[137] -= (fwd_rates[2258] - rev_rates[1116]);
  //sp 149
  sp_rates[148] += (fwd_rates[2258] - rev_rates[1116]);

  //rxn 2259
  //sp 117
  sp_rates[116] += (fwd_rates[2259] - rev_rates[1117]);
  //sp 138
  sp_rates[137] -= (fwd_rates[2259] - rev_rates[1117]);
  //sp 12
  sp_rates[11] += (fwd_rates[2259] - rev_rates[1117]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2259] - rev_rates[1117]);

  //rxn 2260
  //sp 138
  sp_rates[137] -= (fwd_rates[2260] - rev_rates[1118]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2260] - rev_rates[1118]);
  //sp 133
  sp_rates[132] += (fwd_rates[2260] - rev_rates[1118]);
  //sp 149
  sp_rates[148] += (fwd_rates[2260] - rev_rates[1118]);

  //rxn 2261
  //sp 138
  sp_rates[137] -= (fwd_rates[2261] - rev_rates[1119]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2261] - rev_rates[1119]);
  //sp 149
  sp_rates[148] += (fwd_rates[2261] - rev_rates[1119]);
  //sp 136
  sp_rates[135] += (fwd_rates[2261] - rev_rates[1119]);

  //rxn 2262
  //sp 114
  sp_rates[113] -= (fwd_rates[2262] - rev_rates[1120]);
  //sp 12
  sp_rates[11] += (fwd_rates[2262] - rev_rates[1120]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2262] - rev_rates[1120]);
  //sp 1
  (*dy_N) += (fwd_rates[2262] - rev_rates[1120]);

  //rxn 2263
  //sp 114
  sp_rates[113] -= (fwd_rates[2263] - rev_rates[1121]);
  //sp 11
  sp_rates[10] += (fwd_rates[2263] - rev_rates[1121]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2263] - rev_rates[1121]);
  //sp 115
  sp_rates[114] += (fwd_rates[2263] - rev_rates[1121]);

  //rxn 2264
  //sp 12
  sp_rates[11] += (fwd_rates[2264] - rev_rates[1122]);
  //sp 11
  sp_rates[10] -= (fwd_rates[2264] - rev_rates[1122]);
  //sp 1
  (*dy_N) += (fwd_rates[2264] - rev_rates[1122]);
  //sp 115
  sp_rates[114] -= (fwd_rates[2264] - rev_rates[1122]);

  //rxn 2265
  //sp 4
  sp_rates[3] += (fwd_rates[2265] - rev_rates[1123]);
  //sp 151
  sp_rates[150] = -(fwd_rates[2265] - rev_rates[1123]);
  //sp 136
  sp_rates[135] += (fwd_rates[2265] - rev_rates[1123]);

  //rxn 2266
  //sp 3
  sp_rates[2] += (fwd_rates[2266] - rev_rates[1124]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2266] - rev_rates[1124]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2266] - rev_rates[1124]);
  //sp 136
  sp_rates[135] += (fwd_rates[2266] - rev_rates[1124]);

  //rxn 2267
  //sp 136
  sp_rates[135] += (fwd_rates[2267] - rev_rates[1125]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2267] - rev_rates[1125]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2267] - rev_rates[1125]);
  //sp 8
  sp_rates[7] += (fwd_rates[2267] - rev_rates[1125]);

  //rxn 2268
  //sp 136
  sp_rates[135] += (fwd_rates[2268] - rev_rates[1126]);
  //sp 7
  sp_rates[6] += (fwd_rates[2268] - rev_rates[1126]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2268] - rev_rates[1126]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2268] - rev_rates[1126]);

  //rxn 2269
  //sp 136
  sp_rates[135] += (fwd_rates[2269] - rev_rates[1127]);
  //sp 7
  sp_rates[6] += (fwd_rates[2269] - rev_rates[1127]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2269] - rev_rates[1127]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2269] - rev_rates[1127]);

  //rxn 2270
  //sp 25
  sp_rates[24] += (fwd_rates[2270] - rev_rates[1128]);
  //sp 114
  sp_rates[113] += (fwd_rates[2270] - rev_rates[1128]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2270] - rev_rates[1128]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2270] - rev_rates[1128]);

  //rxn 2271
  //sp 127
  sp_rates[126] += (fwd_rates[2271] - rev_rates[1129]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2271] - rev_rates[1129]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2271] - rev_rates[1129]);
  //sp 136
  sp_rates[135] += (fwd_rates[2271] - rev_rates[1129]);

  //rxn 2272
  //sp 131
  sp_rates[130] -= (fwd_rates[2272] - rev_rates[1130]);
  //sp 133
  sp_rates[132] += (fwd_rates[2272] - rev_rates[1130]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2272] - rev_rates[1130]);
  //sp 136
  sp_rates[135] += (fwd_rates[2272] - rev_rates[1130]);

  //rxn 2273
  //sp 129
  sp_rates[128] -= (fwd_rates[2273] - rev_rates[1131]);
  //sp 1
  (*dy_N) += (fwd_rates[2273] - rev_rates[1131]);
  //sp 151
  sp_rates[150] -= (fwd_rates[2273] - rev_rates[1131]);
  //sp 16
  sp_rates[15] += (fwd_rates[2273] - rev_rates[1131]);

  //rxn 2274
  //sp 136
  sp_rates[135] += (fwd_rates[2274] - rev_rates[1132]);
  //sp 4
  sp_rates[3] += (fwd_rates[2274] - rev_rates[1132]);
  //sp 152
  sp_rates[151] = -(fwd_rates[2274] - rev_rates[1132]);

  //rxn 2275
  //sp 151
  sp_rates[150] += (fwd_rates[2275] - rev_rates[1133]);
  //sp 152
  sp_rates[151] -= (fwd_rates[2275] - rev_rates[1133]);

  //rxn 2276
  //sp 136
  sp_rates[135] += (fwd_rates[2276] - rev_rates[1134]);
  //sp 3
  sp_rates[2] += (fwd_rates[2276] - rev_rates[1134]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2276] - rev_rates[1134]);
  //sp 152
  sp_rates[151] -= (fwd_rates[2276] - rev_rates[1134]);

  //rxn 2277
  //sp 138
  sp_rates[137] += (fwd_rates[2277] - rev_rates[1135]);
  //sp 4
  sp_rates[3] += (fwd_rates[2277] - rev_rates[1135]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2277] - rev_rates[1135]);
  //sp 152
  sp_rates[151] -= (fwd_rates[2277] - rev_rates[1135]);

  //rxn 2278
  //sp 136
  sp_rates[135] += (fwd_rates[2278] - rev_rates[1136]);
  //sp 8
  sp_rates[7] += (fwd_rates[2278] - rev_rates[1136]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2278] - rev_rates[1136]);
  //sp 152
  sp_rates[151] -= (fwd_rates[2278] - rev_rates[1136]);

  //rxn 2279
  //sp 152
  sp_rates[151] -= (fwd_rates[2279] - rev_rates[1137]);
  //sp 136
  sp_rates[135] += (fwd_rates[2279] - rev_rates[1137]);
  //sp 7
  sp_rates[6] += (fwd_rates[2279] - rev_rates[1137]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2279] - rev_rates[1137]);

  //rxn 2280
  //sp 136
  sp_rates[135] += (fwd_rates[2280] - rev_rates[1138]);
  //sp 14
  sp_rates[13] += (fwd_rates[2280] - rev_rates[1138]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2280] - rev_rates[1138]);
  //sp 152
  sp_rates[151] -= (fwd_rates[2280] - rev_rates[1138]);

  //rxn 2281
  //sp 129
  sp_rates[128] -= (fwd_rates[2281] - rev_rates[1139]);
  //sp 114
  sp_rates[113] += (fwd_rates[2281] - rev_rates[1139]);
  //sp 11
  sp_rates[10] += (fwd_rates[2281] - rev_rates[1139]);
  //sp 12
  sp_rates[11] -= (fwd_rates[2281] - rev_rates[1139]);

  //rxn 2282
  //sp 11
  sp_rates[10] += (fwd_rates[2282] - rev_rates[1140]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2282] - rev_rates[1140]);
  //sp 12
  sp_rates[11] -= (fwd_rates[2282] - rev_rates[1140]);
  //sp 149
  sp_rates[148] += (fwd_rates[2282] - rev_rates[1140]);

  //rxn 2283
  //sp 129
  sp_rates[128] -= (fwd_rates[2283] - rev_rates[1141]);
  //sp 4
  sp_rates[3] += (fwd_rates[2283] - rev_rates[1141]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2283] - rev_rates[1141]);
  //sp 152
  sp_rates[151] += (fwd_rates[2283] - rev_rates[1141]);

  //rxn 2284
  //sp 129
  sp_rates[128] -= (fwd_rates[2284] - rev_rates[1142]);
  //sp 151
  sp_rates[150] += (fwd_rates[2284] - rev_rates[1142]);
  //sp 4
  sp_rates[3] += (fwd_rates[2284] - rev_rates[1142]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2284] - rev_rates[1142]);

  //rxn 2285
  //sp 136
  sp_rates[135] += (fwd_rates[2285] - rev_rates[1143]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2285] - rev_rates[1143]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2285] - rev_rates[1143]);
  //sp 7
  sp_rates[6] += (fwd_rates[2285] - rev_rates[1143]);

  //rxn 2286
  //sp 114
  sp_rates[113] -= (fwd_rates[2286] - rev_rates[1144]);
  //sp 151
  sp_rates[150] += (fwd_rates[2286] - rev_rates[1144]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2286] - rev_rates[1144]);
  //sp 8
  sp_rates[7] += (fwd_rates[2286] - rev_rates[1144]);

  //rxn 2287
  //sp 156
  sp_rates[155] = (fwd_rates[2287] - rev_rates[1145]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2287] - rev_rates[1145]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2287] - rev_rates[1145]);

  //rxn 2288
  //sp 157
  sp_rates[156] = (fwd_rates[2288] - rev_rates[1146]);
  //sp 4
  sp_rates[3] += (fwd_rates[2288] - rev_rates[1146]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2288] - rev_rates[1146]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2288] - rev_rates[1146]);

  //rxn 2289
  //sp 4
  sp_rates[3] += (fwd_rates[2289] - rev_rates[1147]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2289] - rev_rates[1147]);
  //sp 158
  sp_rates[157] = (fwd_rates[2289] - rev_rates[1147]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2289] - rev_rates[1147]);

  //rxn 2290
  //sp 3
  sp_rates[2] += (fwd_rates[2290] - rev_rates[1148]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2290] - rev_rates[1148]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2290] - rev_rates[1148]);
  //sp 155
  sp_rates[154] = (fwd_rates[2290] - rev_rates[1148]);

  //rxn 2291
  //sp 131
  sp_rates[130] -= (fwd_rates[2291] - rev_rates[1149]);
  //sp 4
  sp_rates[3] += (fwd_rates[2291] - rev_rates[1149]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2291] - rev_rates[1149]);
  //sp 155
  sp_rates[154] += (fwd_rates[2291] - rev_rates[1149]);

  //rxn 2292
  //sp 129
  sp_rates[128] += (fwd_rates[2292] - rev_rates[1150]);
  //sp 131
  sp_rates[130] -= (fwd_rates[2292] - rev_rates[1150]);
  //sp 14
  sp_rates[13] += (fwd_rates[2292] - rev_rates[1150]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2292] - rev_rates[1150]);

  //rxn 2293
  //sp 3
  sp_rates[2] += (fwd_rates[2293] - rev_rates[1151]);
  //sp 155
  sp_rates[154] -= (fwd_rates[2293] - rev_rates[1151]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2293] - rev_rates[1151]);
  //sp 151
  sp_rates[150] += (fwd_rates[2293] - rev_rates[1151]);

  //rxn 2294
  //sp 155
  sp_rates[154] -= (fwd_rates[2294] - rev_rates[1152]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2294] - rev_rates[1152]);
  //sp 151
  sp_rates[150] += (fwd_rates[2294] - rev_rates[1152]);
  //sp 8
  sp_rates[7] += (fwd_rates[2294] - rev_rates[1152]);

  //rxn 2295
  //sp 155
  sp_rates[154] -= (fwd_rates[2295] - rev_rates[1153]);
  //sp 151
  sp_rates[150] += (fwd_rates[2295] - rev_rates[1153]);
  //sp 7
  sp_rates[6] += (fwd_rates[2295] - rev_rates[1153]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2295] - rev_rates[1153]);

  //rxn 2296
  //sp 155
  sp_rates[154] -= (fwd_rates[2296] - rev_rates[1154]);
  //sp 151
  sp_rates[150] += (fwd_rates[2296] - rev_rates[1154]);
  //sp 14
  sp_rates[13] += (fwd_rates[2296] - rev_rates[1154]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2296] - rev_rates[1154]);

  //rxn 2297
  //sp 155
  sp_rates[154] -= (fwd_rates[2297] - rev_rates[1155]);
  //sp 151
  sp_rates[150] += (fwd_rates[2297] - rev_rates[1155]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2297] - rev_rates[1155]);
  //sp 127
  sp_rates[126] += (fwd_rates[2297] - rev_rates[1155]);

  //rxn 2298
  //sp 3
  sp_rates[2] += (fwd_rates[2298] - rev_rates[1156]);
  //sp 155
  sp_rates[154] -= (fwd_rates[2298] - rev_rates[1156]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2298] - rev_rates[1156]);
  //sp 152
  sp_rates[151] += (fwd_rates[2298] - rev_rates[1156]);

  //rxn 2299
  //sp 152
  sp_rates[151] += (fwd_rates[2299] - rev_rates[1157]);
  //sp 155
  sp_rates[154] -= (fwd_rates[2299] - rev_rates[1157]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2299] - rev_rates[1157]);
  //sp 8
  sp_rates[7] += (fwd_rates[2299] - rev_rates[1157]);

  //rxn 2300
  //sp 152
  sp_rates[151] += (fwd_rates[2300] - rev_rates[1158]);
  //sp 155
  sp_rates[154] -= (fwd_rates[2300] - rev_rates[1158]);
  //sp 7
  sp_rates[6] += (fwd_rates[2300] - rev_rates[1158]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2300] - rev_rates[1158]);

  //rxn 2301
  //sp 155
  sp_rates[154] -= (fwd_rates[2301] - rev_rates[1159]);
  //sp 14
  sp_rates[13] += (fwd_rates[2301] - rev_rates[1159]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2301] - rev_rates[1159]);
  //sp 152
  sp_rates[151] += (fwd_rates[2301] - rev_rates[1159]);

  //rxn 2302
  //sp 155
  sp_rates[154] -= (fwd_rates[2302] - rev_rates[1160]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2302] - rev_rates[1160]);
  //sp 127
  sp_rates[126] += (fwd_rates[2302] - rev_rates[1160]);
  //sp 152
  sp_rates[151] += (fwd_rates[2302] - rev_rates[1160]);

  //rxn 2303
  //sp 25
  sp_rates[24] += (fwd_rates[2303] - rev_rates[1161]);
  //sp 155
  sp_rates[154] -= (fwd_rates[2303] - rev_rates[1161]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2303] - rev_rates[1161]);
  //sp 131
  sp_rates[130] += (fwd_rates[2303] - rev_rates[1161]);

  //rxn 2304
  //sp 155
  sp_rates[154] += (fwd_rates[2304] - rev_rates[1162]);
  //sp 4
  sp_rates[3] += (fwd_rates[2304] - rev_rates[1162]);
  //sp 158
  sp_rates[157] -= (fwd_rates[2304] - rev_rates[1162]);

  //rxn 2305
  //sp 3
  sp_rates[2] += (fwd_rates[2305] - rev_rates[1163]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2305] - rev_rates[1163]);
  //sp 158
  sp_rates[157] -= (fwd_rates[2305] - rev_rates[1163]);
  //sp 155
  sp_rates[154] += (fwd_rates[2305] - rev_rates[1163]);

  //rxn 2306
  //sp 155
  sp_rates[154] += (fwd_rates[2306] - rev_rates[1164]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2306] - rev_rates[1164]);
  //sp 158
  sp_rates[157] -= (fwd_rates[2306] - rev_rates[1164]);
  //sp 8
  sp_rates[7] += (fwd_rates[2306] - rev_rates[1164]);

  //rxn 2307
  //sp 155
  sp_rates[154] += (fwd_rates[2307] - rev_rates[1165]);
  //sp 158
  sp_rates[157] -= (fwd_rates[2307] - rev_rates[1165]);
  //sp 7
  sp_rates[6] += (fwd_rates[2307] - rev_rates[1165]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2307] - rev_rates[1165]);

  //rxn 2308
  //sp 14
  sp_rates[13] += (fwd_rates[2308] - rev_rates[1166]);
  //sp 158
  sp_rates[157] -= (fwd_rates[2308] - rev_rates[1166]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2308] - rev_rates[1166]);
  //sp 155
  sp_rates[154] += (fwd_rates[2308] - rev_rates[1166]);

  //rxn 2309
  //sp 155
  sp_rates[154] += (fwd_rates[2309] - rev_rates[1167]);
  //sp 4
  sp_rates[3] += (fwd_rates[2309] - rev_rates[1167]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2309] - rev_rates[1167]);

  //rxn 2310
  //sp 157
  sp_rates[156] -= (fwd_rates[2310] - rev_rates[1168]);
  //sp 155
  sp_rates[154] += (fwd_rates[2310] - rev_rates[1168]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2310] - rev_rates[1168]);
  //sp 10
  sp_rates[9] += (fwd_rates[2310] - rev_rates[1168]);

  //rxn 2311
  //sp 5
  sp_rates[4] -= (fwd_rates[2311] - rev_rates[1169]);
  //sp 6
  sp_rates[5] += (fwd_rates[2311] - rev_rates[1169]);
  //sp 133
  sp_rates[132] += (fwd_rates[2311] - rev_rates[1169]);
  //sp 25
  sp_rates[24] += (fwd_rates[2311] - rev_rates[1169]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2311] - rev_rates[1169]);

  //rxn 2312
  //sp 3
  sp_rates[2] += (fwd_rates[2312] - rev_rates[1170]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2312] - rev_rates[1170]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2312] - rev_rates[1170]);
  //sp 155
  sp_rates[154] += (fwd_rates[2312] - rev_rates[1170]);

  //rxn 2313
  //sp 25
  sp_rates[24] += (fwd_rates[2313] - rev_rates[1171]);
  //sp 133
  sp_rates[132] += (fwd_rates[2313] - rev_rates[1171]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2313] - rev_rates[1171]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2313] - rev_rates[1171]);

  //rxn 2314
  //sp 155
  sp_rates[154] += (fwd_rates[2314] - rev_rates[1172]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2314] - rev_rates[1172]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2314] - rev_rates[1172]);
  //sp 8
  sp_rates[7] += (fwd_rates[2314] - rev_rates[1172]);

  //rxn 2315
  //sp 133
  sp_rates[132] += (fwd_rates[2315] - rev_rates[1173]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2315] - rev_rates[1173]);
  //sp 24
  sp_rates[23] += (fwd_rates[2315] - rev_rates[1173]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2315] - rev_rates[1173]);

  //rxn 2316
  //sp 155
  sp_rates[154] += (fwd_rates[2316] - rev_rates[1174]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2316] - rev_rates[1174]);
  //sp 7
  sp_rates[6] += (fwd_rates[2316] - rev_rates[1174]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2316] - rev_rates[1174]);

  //rxn 2317
  //sp 133
  sp_rates[132] += (fwd_rates[2317] - rev_rates[1175]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2317] - rev_rates[1175]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2317] - rev_rates[1175]);
  //sp 31
  sp_rates[30] += (fwd_rates[2317] - rev_rates[1175]);

  //rxn 2318
  //sp 155
  sp_rates[154] += (fwd_rates[2318] - rev_rates[1176]);
  //sp 157
  sp_rates[156] -= (fwd_rates[2318] - rev_rates[1176]);
  //sp 14
  sp_rates[13] += (fwd_rates[2318] - rev_rates[1176]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2318] - rev_rates[1176]);

  //rxn 2319
  //sp 157
  sp_rates[156] += (fwd_rates[2319] - rev_rates[1177]);
  //sp 3
  sp_rates[2] += (fwd_rates[2319] - rev_rates[1177]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2319] - rev_rates[1177]);
  //sp 156
  sp_rates[155] -= (fwd_rates[2319] - rev_rates[1177]);

  //rxn 2320
  //sp 156
  sp_rates[155] -= (fwd_rates[2320] - rev_rates[1178]);
  //sp 157
  sp_rates[156] += (fwd_rates[2320] - rev_rates[1178]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2320] - rev_rates[1178]);
  //sp 8
  sp_rates[7] += (fwd_rates[2320] - rev_rates[1178]);

  //rxn 2321
  //sp 156
  sp_rates[155] -= (fwd_rates[2321] - rev_rates[1179]);
  //sp 157
  sp_rates[156] += (fwd_rates[2321] - rev_rates[1179]);
  //sp 7
  sp_rates[6] += (fwd_rates[2321] - rev_rates[1179]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2321] - rev_rates[1179]);

  //rxn 2322
  //sp 156
  sp_rates[155] -= (fwd_rates[2322] - rev_rates[1180]);
  //sp 157
  sp_rates[156] += (fwd_rates[2322] - rev_rates[1180]);
  //sp 14
  sp_rates[13] += (fwd_rates[2322] - rev_rates[1180]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2322] - rev_rates[1180]);

  //rxn 2323
  //sp 157
  sp_rates[156] += (fwd_rates[2323] - rev_rates[1181]);
  //sp 156
  sp_rates[155] -= (fwd_rates[2323] - rev_rates[1181]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2323] - rev_rates[1181]);
  //sp 127
  sp_rates[126] += (fwd_rates[2323] - rev_rates[1181]);

  //rxn 2324
  //sp 3
  sp_rates[2] += (fwd_rates[2324] - rev_rates[1182]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2324] - rev_rates[1182]);
  //sp 156
  sp_rates[155] -= (fwd_rates[2324] - rev_rates[1182]);
  //sp 158
  sp_rates[157] += (fwd_rates[2324] - rev_rates[1182]);

  //rxn 2325
  //sp 156
  sp_rates[155] -= (fwd_rates[2325] - rev_rates[1183]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2325] - rev_rates[1183]);
  //sp 158
  sp_rates[157] += (fwd_rates[2325] - rev_rates[1183]);
  //sp 8
  sp_rates[7] += (fwd_rates[2325] - rev_rates[1183]);

  //rxn 2326
  //sp 156
  sp_rates[155] -= (fwd_rates[2326] - rev_rates[1184]);
  //sp 158
  sp_rates[157] += (fwd_rates[2326] - rev_rates[1184]);
  //sp 7
  sp_rates[6] += (fwd_rates[2326] - rev_rates[1184]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2326] - rev_rates[1184]);

  //rxn 2327
  //sp 156
  sp_rates[155] -= (fwd_rates[2327] - rev_rates[1185]);
  //sp 14
  sp_rates[13] += (fwd_rates[2327] - rev_rates[1185]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2327] - rev_rates[1185]);
  //sp 158
  sp_rates[157] += (fwd_rates[2327] - rev_rates[1185]);

  //rxn 2328
  //sp 156
  sp_rates[155] -= (fwd_rates[2328] - rev_rates[1186]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2328] - rev_rates[1186]);
  //sp 158
  sp_rates[157] += (fwd_rates[2328] - rev_rates[1186]);
  //sp 127
  sp_rates[126] += (fwd_rates[2328] - rev_rates[1186]);

  //rxn 2329
  //sp 131
  sp_rates[130] += (fwd_rates[2329] - rev_rates[1187]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2329] - rev_rates[1187]);
  //sp 14
  sp_rates[13] += (fwd_rates[2329] - rev_rates[1187]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2329] - rev_rates[1187]);

  //rxn 2330
  //sp 127
  sp_rates[126] += (fwd_rates[2330] - rev_rates[1188]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2330] - rev_rates[1188]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2330] - rev_rates[1188]);
  //sp 16
  sp_rates[15] += (fwd_rates[2330] - rev_rates[1188]);

  //rxn 2331
  //sp 129
  sp_rates[128] -= (fwd_rates[2331] - rev_rates[1189]);
  //sp 136
  sp_rates[135] += (fwd_rates[2331] - rev_rates[1189]);
  //sp 4
  sp_rates[3] += (fwd_rates[2331] - rev_rates[1189]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2331] - rev_rates[1189]);

  //rxn 2332
  //sp 136
  sp_rates[135] += (fwd_rates[2332] - rev_rates[1190]);
  //sp 131
  sp_rates[130] += (fwd_rates[2332] - rev_rates[1190]);
  //sp 1
  (*dy_N) -= (fwd_rates[2332] - rev_rates[1190]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2332] - rev_rates[1190]);

  //rxn 2333
  //sp 129
  sp_rates[128] -= (fwd_rates[2333] - rev_rates[1191]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2333] - rev_rates[1191]);
  //sp 4
  sp_rates[3] += (fwd_rates[2333] - rev_rates[1191]);
  //sp 147
  sp_rates[146] += (fwd_rates[2333] - rev_rates[1191]);

  //rxn 2334
  //sp 114
  sp_rates[113] += (fwd_rates[2334] - rev_rates[1192]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2334] - rev_rates[1192]);
  //sp 115
  sp_rates[114] -= (fwd_rates[2334] - rev_rates[1192]);
  //sp 136
  sp_rates[135] += (fwd_rates[2334] - rev_rates[1192]);

  //rxn 2335
  //sp 129
  sp_rates[128] += (fwd_rates[2335] - rev_rates[1193]);
  //sp 18
  sp_rates[17] -= (fwd_rates[2335] - rev_rates[1193]);
  //sp 147
  sp_rates[146] += (fwd_rates[2335] - rev_rates[1193]);
  //sp 1
  (*dy_N) -= (fwd_rates[2335] - rev_rates[1193]);

  //rxn 2336
  //sp 18
  sp_rates[17] -= (fwd_rates[2336] - rev_rates[1194]);
  //sp 147
  sp_rates[146] += (fwd_rates[2336] - rev_rates[1194]);
  //sp 6
  sp_rates[5] += (fwd_rates[2336] - rev_rates[1194]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2336] - rev_rates[1194]);

  //rxn 2337
  //sp 129
  sp_rates[128] += (fwd_rates[2337] - rev_rates[1195]);
  //sp 18
  sp_rates[17] -= (fwd_rates[2337] - rev_rates[1195]);
  //sp 11
  sp_rates[10] += (fwd_rates[2337] - rev_rates[1195]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2337] - rev_rates[1195]);

  //rxn 2338
  //sp 114
  sp_rates[113] -= (fwd_rates[2338] - rev_rates[1196]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2338] - rev_rates[1196]);
  //sp 6
  sp_rates[5] += (fwd_rates[2338] - rev_rates[1196]);
  //sp 136
  sp_rates[135] += (fwd_rates[2338] - rev_rates[1196]);

  //rxn 2339
  //sp 114
  sp_rates[113] -= (fwd_rates[2339] - rev_rates[1197]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2339] - rev_rates[1197]);
  //sp 4
  sp_rates[3] += (fwd_rates[2339] - rev_rates[1197]);
  //sp 149
  sp_rates[148] += (fwd_rates[2339] - rev_rates[1197]);

  //rxn 2340
  //sp 129
  sp_rates[128] += (fwd_rates[2340] - rev_rates[1198]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2340] - rev_rates[1198]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2340] - rev_rates[1198]);
  //sp 26
  sp_rates[25] += (fwd_rates[2340] - rev_rates[1198]);

  //rxn 2341
  //sp 114
  sp_rates[113] -= (fwd_rates[2341] - rev_rates[1199]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2341] - rev_rates[1199]);
  //sp 11
  sp_rates[10] += (fwd_rates[2341] - rev_rates[1199]);
  //sp 131
  sp_rates[130] += (fwd_rates[2341] - rev_rates[1199]);

  //rxn 2342
  //sp 114
  sp_rates[113] -= (fwd_rates[2342] - rev_rates[1200]);
  //sp 19
  sp_rates[18] -= (fwd_rates[2342] - rev_rates[1200]);
  //sp 147
  sp_rates[146] += (fwd_rates[2342] - rev_rates[1200]);
  //sp 8
  sp_rates[7] += (fwd_rates[2342] - rev_rates[1200]);

  //rxn 2343
  //sp 136
  sp_rates[135] += (fwd_rates[2343] - rev_rates[1201]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2343] - rev_rates[1201]);
  //sp 8
  sp_rates[7] += (fwd_rates[2343] - rev_rates[1201]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2343] - rev_rates[1201]);

  //rxn 2344
  //sp 114
  sp_rates[113] -= (fwd_rates[2344] - rev_rates[1202]);
  //sp 138
  sp_rates[137] += (fwd_rates[2344] - rev_rates[1202]);
  //sp 4
  sp_rates[3] += (fwd_rates[2344] - rev_rates[1202]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2344] - rev_rates[1202]);

  //rxn 2345
  //sp 114
  sp_rates[113] -= (fwd_rates[2345] - rev_rates[1203]);
  //sp 139
  sp_rates[138] += (fwd_rates[2345] - rev_rates[1203]);
  //sp 4
  sp_rates[3] += (fwd_rates[2345] - rev_rates[1203]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2345] - rev_rates[1203]);

  //rxn 2346
  //sp 114
  sp_rates[113] -= (fwd_rates[2346] - rev_rates[1204]);
  //sp 11
  sp_rates[10] += (fwd_rates[2346] - rev_rates[1204]);
  //sp 133
  sp_rates[132] += (fwd_rates[2346] - rev_rates[1204]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2346] - rev_rates[1204]);

  //rxn 2347
  //sp 114
  sp_rates[113] -= (fwd_rates[2347] - rev_rates[1205]);
  //sp 6
  sp_rates[5] += (fwd_rates[2347] - rev_rates[1205]);
  //sp 151
  sp_rates[150] += (fwd_rates[2347] - rev_rates[1205]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2347] - rev_rates[1205]);

  //rxn 2348
  //sp 17
  sp_rates[16] -= (fwd_rates[2348] - rev_rates[1206]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2348] - rev_rates[1206]);
  //sp 136
  sp_rates[135] += (fwd_rates[2348] - rev_rates[1206]);
  //sp 8
  sp_rates[7] += (fwd_rates[2348] - rev_rates[1206]);

  //rxn 2349
  //sp 17
  sp_rates[16] -= (fwd_rates[2349] - rev_rates[1207]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2349] - rev_rates[1207]);
  //sp 138
  sp_rates[137] += (fwd_rates[2349] - rev_rates[1207]);
  //sp 4
  sp_rates[3] += (fwd_rates[2349] - rev_rates[1207]);

  //rxn 2350
  //sp 17
  sp_rates[16] -= (fwd_rates[2350] - rev_rates[1208]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2350] - rev_rates[1208]);
  //sp 139
  sp_rates[138] += (fwd_rates[2350] - rev_rates[1208]);
  //sp 4
  sp_rates[3] += (fwd_rates[2350] - rev_rates[1208]);

  //rxn 2351
  //sp 17
  sp_rates[16] -= (fwd_rates[2351] - rev_rates[1209]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2351] - rev_rates[1209]);
  //sp 11
  sp_rates[10] += (fwd_rates[2351] - rev_rates[1209]);
  //sp 133
  sp_rates[132] += (fwd_rates[2351] - rev_rates[1209]);

  //rxn 2352
  //sp 17
  sp_rates[16] -= (fwd_rates[2352] - rev_rates[1210]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2352] - rev_rates[1210]);
  //sp 6
  sp_rates[5] += (fwd_rates[2352] - rev_rates[1210]);
  //sp 151
  sp_rates[150] += (fwd_rates[2352] - rev_rates[1210]);

  //rxn 2353
  //sp 17
  sp_rates[16] -= (fwd_rates[2353] - rev_rates[1211]);
  //sp 16
  sp_rates[15] += (fwd_rates[2353] - rev_rates[1211]);

  //rxn 2354
  //sp 139
  sp_rates[138] -= (fwd_rates[2354] - rev_rates[1212]);
  //sp 6
  sp_rates[5] += (fwd_rates[2354] - rev_rates[1212]);
  //sp 136
  sp_rates[135] += (fwd_rates[2354] - rev_rates[1212]);

  //rxn 2355
  //sp 137
  sp_rates[136] = (fwd_rates[2355] - rev_rates[1213]);
  //sp 136
  sp_rates[135] -= (fwd_rates[2355] - rev_rates[1213]);

  //rxn 2356
  //sp 137
  sp_rates[136] -= (fwd_rates[2356] - rev_rates[1214]);
  //sp 136
  sp_rates[135] += (fwd_rates[2356] - rev_rates[1214]);

  //rxn 2357
  //sp 137
  sp_rates[136] -= (fwd_rates[2357] - rev_rates[1215]);
  //sp 11
  sp_rates[10] += (fwd_rates[2357] - rev_rates[1215]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2357] - rev_rates[1215]);
  //sp 131
  sp_rates[130] += (fwd_rates[2357] - rev_rates[1215]);

  //rxn 2358
  //sp 137
  sp_rates[136] -= (fwd_rates[2358] - rev_rates[1216]);
  //sp 138
  sp_rates[137] += (fwd_rates[2358] - rev_rates[1216]);
  //sp 4
  sp_rates[3] += (fwd_rates[2358] - rev_rates[1216]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2358] - rev_rates[1216]);

  //rxn 2359
  //sp 137
  sp_rates[136] -= (fwd_rates[2359] - rev_rates[1217]);
  //sp 147
  sp_rates[146] += (fwd_rates[2359] - rev_rates[1217]);
  //sp 7
  sp_rates[6] += (fwd_rates[2359] - rev_rates[1217]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2359] - rev_rates[1217]);

  //rxn 2360
  //sp 137
  sp_rates[136] -= (fwd_rates[2360] - rev_rates[1218]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2360] - rev_rates[1218]);
  //sp 4
  sp_rates[3] += (fwd_rates[2360] - rev_rates[1218]);
  //sp 153
  sp_rates[152] = (fwd_rates[2360] - rev_rates[1218]);

  //rxn 2361
  //sp 147
  sp_rates[146] -= (fwd_rates[2361] - rev_rates[1219]);
  //sp 30
  sp_rates[29] -= (fwd_rates[2361] - rev_rates[1219]);
  //sp 31
  sp_rates[30] += (fwd_rates[2361] - rev_rates[1219]);
  //sp 136
  sp_rates[135] += (fwd_rates[2361] - rev_rates[1219]);

  //rxn 2362
  //sp 138
  sp_rates[137] += (fwd_rates[2362] - rev_rates[1220]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2362] - rev_rates[1220]);
  //sp 30
  sp_rates[29] -= (fwd_rates[2362] - rev_rates[1220]);
  //sp 31
  sp_rates[30] += (fwd_rates[2362] - rev_rates[1220]);

  //rxn 2363
  //sp 34
  sp_rates[33] -= (fwd_rates[2363] - rev_rates[1221]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2363] - rev_rates[1221]);
  //sp 35
  sp_rates[34] += (fwd_rates[2363] - rev_rates[1221]);
  //sp 136
  sp_rates[135] += (fwd_rates[2363] - rev_rates[1221]);

  //rxn 2364
  //sp 147
  sp_rates[146] -= (fwd_rates[2364] - rev_rates[1222]);
  //sp 62
  sp_rates[61] -= (fwd_rates[2364] - rev_rates[1222]);
  //sp 63
  sp_rates[62] += (fwd_rates[2364] - rev_rates[1222]);
  //sp 136
  sp_rates[135] += (fwd_rates[2364] - rev_rates[1222]);

  //rxn 2365
  //sp 114
  sp_rates[113] -= (fwd_rates[2365] - rev_rates[1223]);
  //sp 35
  sp_rates[34] -= (fwd_rates[2365] - rev_rates[1223]);
  //sp 36
  sp_rates[35] += (fwd_rates[2365] - rev_rates[1223]);
  //sp 117
  sp_rates[116] += (fwd_rates[2365] - rev_rates[1223]);

  //rxn 2366
  //sp 25
  sp_rates[24] += (fwd_rates[2366] - rev_rates[1224]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2366] - rev_rates[1224]);
  //sp 35
  sp_rates[34] -= (fwd_rates[2366] - rev_rates[1224]);
  //sp 136
  sp_rates[135] += (fwd_rates[2366] - rev_rates[1224]);

  //rxn 2367
  //sp 49
  sp_rates[48] += (fwd_rates[2367] - rev_rates[1225]);
  //sp 36
  sp_rates[35] -= (fwd_rates[2367] - rev_rates[1225]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2367] - rev_rates[1225]);
  //sp 136
  sp_rates[135] += (fwd_rates[2367] - rev_rates[1225]);

  //rxn 2368
  //sp 49
  sp_rates[48] += (fwd_rates[2368] - rev_rates[1226]);
  //sp 136
  sp_rates[135] += (fwd_rates[2368] - rev_rates[1226]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2368] - rev_rates[1226]);
  //sp 48
  sp_rates[47] -= (fwd_rates[2368] - rev_rates[1226]);

  //rxn 2369
  //sp 147
  sp_rates[146] -= (fwd_rates[2369] - rev_rates[1227]);
  //sp 14
  sp_rates[13] -= (fwd_rates[2369] - rev_rates[1227]);
  //sp 15
  sp_rates[14] += (fwd_rates[2369] - rev_rates[1227]);
  //sp 136
  sp_rates[135] += (fwd_rates[2369] - rev_rates[1227]);

  //rxn 2370
  //sp 146
  sp_rates[145] = -(fwd_rates[2370] - rev_rates[1228]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2370] - rev_rates[1228]);
  //sp 15
  sp_rates[14] += (fwd_rates[2370] - rev_rates[1228]);
  //sp 136
  sp_rates[135] += (fwd_rates[2370] - rev_rates[1228]);

  //rxn 2371
  //sp 137
  sp_rates[136] += (fwd_rates[2371] - rev_rates[1229]);
  //sp 146
  sp_rates[145] -= (fwd_rates[2371] - rev_rates[1229]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2371] - rev_rates[1229]);
  //sp 15
  sp_rates[14] += (fwd_rates[2371] - rev_rates[1229]);

  //rxn 2372
  //sp 154
  sp_rates[153] = (fwd_rates[2372] - rev_rates[1230]);
  //sp 146
  sp_rates[145] -= (fwd_rates[2372] - rev_rates[1230]);
  //sp 3
  sp_rates[2] += (fwd_rates[2372] - rev_rates[1230]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2372] - rev_rates[1230]);

  //rxn 2373
  //sp 146
  sp_rates[145] -= (fwd_rates[2373] - rev_rates[1231]);
  //sp 35
  sp_rates[34] -= (fwd_rates[2373] - rev_rates[1231]);
  //sp 154
  sp_rates[153] += (fwd_rates[2373] - rev_rates[1231]);
  //sp 34
  sp_rates[33] += (fwd_rates[2373] - rev_rates[1231]);

  //rxn 2374
  //sp 146
  sp_rates[145] -= (fwd_rates[2374] - rev_rates[1232]);
  //sp 149
  sp_rates[148] += (fwd_rates[2374] - rev_rates[1232]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2374] - rev_rates[1232]);
  //sp 15
  sp_rates[14] += (fwd_rates[2374] - rev_rates[1232]);

  //rxn 2375
  //sp 146
  sp_rates[145] -= (fwd_rates[2375] - rev_rates[1233]);
  //sp 154
  sp_rates[153] += (fwd_rates[2375] - rev_rates[1233]);
  //sp 7
  sp_rates[6] += (fwd_rates[2375] - rev_rates[1233]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2375] - rev_rates[1233]);

  //rxn 2376
  //sp 146
  sp_rates[145] -= (fwd_rates[2376] - rev_rates[1234]);
  //sp 154
  sp_rates[153] += (fwd_rates[2376] - rev_rates[1234]);
  //sp 4
  sp_rates[3] += (fwd_rates[2376] - rev_rates[1234]);

  //rxn 2377
  //sp 146
  sp_rates[145] -= (fwd_rates[2377] - rev_rates[1235]);
  //sp 147
  sp_rates[146] += (fwd_rates[2377] - rev_rates[1235]);
  //sp 15
  sp_rates[14] += (fwd_rates[2377] - rev_rates[1235]);

  //rxn 2378
  //sp 146
  sp_rates[145] += (fwd_rates[2378] - rev_rates[1236]);
  //sp 140
  sp_rates[139] -= (fwd_rates[2378] - rev_rates[1236]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2378] - rev_rates[1236]);
  //sp 8
  sp_rates[7] += (fwd_rates[2378] - rev_rates[1236]);

  //rxn 2379
  //sp 25
  sp_rates[24] += (fwd_rates[2379] - rev_rates[1237]);
  //sp 154
  sp_rates[153] -= (fwd_rates[2379] - rev_rates[1237]);
  //sp 147
  sp_rates[146] += (fwd_rates[2379] - rev_rates[1237]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2379] - rev_rates[1237]);

  //rxn 2380
  //sp 154
  sp_rates[153] += (fwd_rates[2380] - rev_rates[1238]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2380] - rev_rates[1238]);
  //sp 8
  sp_rates[7] += (fwd_rates[2380] - rev_rates[1238]);
  //sp 24
  sp_rates[23] -= (fwd_rates[2380] - rev_rates[1238]);

  //rxn 2381
  //sp 154
  sp_rates[153] += (fwd_rates[2381] - rev_rates[1239]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2381] - rev_rates[1239]);
  //sp 4
  sp_rates[3] += (fwd_rates[2381] - rev_rates[1239]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2381] - rev_rates[1239]);

  //rxn 2382
  //sp 146
  sp_rates[145] -= (fwd_rates[2382] - rev_rates[1240]);
  //sp 154
  sp_rates[153] += (fwd_rates[2382] - rev_rates[1240]);
  //sp 14
  sp_rates[13] += (fwd_rates[2382] - rev_rates[1240]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2382] - rev_rates[1240]);

  //rxn 2383
  //sp 146
  sp_rates[145] -= (fwd_rates[2383] - rev_rates[1241]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2383] - rev_rates[1241]);
  //sp 154
  sp_rates[153] += (fwd_rates[2383] - rev_rates[1241]);
  //sp 136
  sp_rates[135] += (fwd_rates[2383] - rev_rates[1241]);

  //rxn 2384
  //sp 153
  sp_rates[152] += (fwd_rates[2384] - rev_rates[1242]);
  //sp 146
  sp_rates[145] -= (fwd_rates[2384] - rev_rates[1242]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2384] - rev_rates[1242]);
  //sp 15
  sp_rates[14] += (fwd_rates[2384] - rev_rates[1242]);

  //rxn 2385
  //sp 129
  sp_rates[128] -= (fwd_rates[2385] - rev_rates[1243]);
  //sp 136
  sp_rates[135] += (fwd_rates[2385] - rev_rates[1243]);
  //sp 35
  sp_rates[34] -= (fwd_rates[2385] - rev_rates[1243]);
  //sp 16
  sp_rates[15] += (fwd_rates[2385] - rev_rates[1243]);

  //rxn 2386
  //sp 49
  sp_rates[48] -= (fwd_rates[2386] - rev_rates[1244]);
  //sp 11
  sp_rates[10] += (fwd_rates[2386] - rev_rates[1244]);
  //sp 129
  sp_rates[128] -= (fwd_rates[2386] - rev_rates[1244]);
  //sp 136
  sp_rates[135] += (fwd_rates[2386] - rev_rates[1244]);

  //rxn 2387
  //sp 49
  sp_rates[48] -= (fwd_rates[2387] - rev_rates[1245]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2387] - rev_rates[1245]);
  //sp 11
  sp_rates[10] += (fwd_rates[2387] - rev_rates[1245]);
  //sp 139
  sp_rates[138] += (fwd_rates[2387] - rev_rates[1245]);

  //rxn 2388
  //sp 49
  sp_rates[48] -= (fwd_rates[2388] - rev_rates[1246]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2388] - rev_rates[1246]);
  //sp 12
  sp_rates[11] += (fwd_rates[2388] - rev_rates[1246]);
  //sp 136
  sp_rates[135] += (fwd_rates[2388] - rev_rates[1246]);

  //rxn 2389
  //sp 49
  sp_rates[48] -= (fwd_rates[2389] - rev_rates[1247]);
  //sp 12
  sp_rates[11] += (fwd_rates[2389] - rev_rates[1247]);
  //sp 139
  sp_rates[138] += (fwd_rates[2389] - rev_rates[1247]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2389] - rev_rates[1247]);

  //rxn 2390
  //sp 114
  sp_rates[113] -= (fwd_rates[2390] - rev_rates[1248]);
  //sp 11
  sp_rates[10] += (fwd_rates[2390] - rev_rates[1248]);
  //sp 37
  sp_rates[36] -= (fwd_rates[2390] - rev_rates[1248]);
  //sp 136
  sp_rates[135] += (fwd_rates[2390] - rev_rates[1248]);

  //rxn 2391
  //sp 114
  sp_rates[113] -= (fwd_rates[2391] - rev_rates[1249]);
  //sp 147
  sp_rates[146] += (fwd_rates[2391] - rev_rates[1249]);
  //sp 37
  sp_rates[36] -= (fwd_rates[2391] - rev_rates[1249]);
  //sp 26
  sp_rates[25] += (fwd_rates[2391] - rev_rates[1249]);

  //rxn 2392
  //sp 129
  sp_rates[128] -= (fwd_rates[2392] - rev_rates[1250]);
  //sp 136
  sp_rates[135] += (fwd_rates[2392] - rev_rates[1250]);
  //sp 36
  sp_rates[35] += (fwd_rates[2392] - rev_rates[1250]);
  //sp 80
  sp_rates[79] -= (fwd_rates[2392] - rev_rates[1250]);

  //rxn 2393
  //sp 115
  sp_rates[114] -= (fwd_rates[2393] - rev_rates[1251]) * pres_mod[59];
  //sp 6
  sp_rates[5] += (fwd_rates[2393] - rev_rates[1251]) * pres_mod[59];
  //sp 1
  (*dy_N) += (fwd_rates[2393] - rev_rates[1251]) * pres_mod[59];

  //rxn 2394
  //sp 115
  sp_rates[114] -= (fwd_rates[2394] - rev_rates[1252]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2394] - rev_rates[1252]);
  //sp 1
  (*dy_N) += (fwd_rates[2394] - rev_rates[1252]);
  //sp 8
  sp_rates[7] += (fwd_rates[2394] - rev_rates[1252]);

  //rxn 2395
  //sp 115
  sp_rates[114] -= (fwd_rates[2395] - rev_rates[1253]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2395] - rev_rates[1253]);
  //sp 1
  (*dy_N) += (fwd_rates[2395] - rev_rates[1253]);
  //sp 8
  sp_rates[7] += (fwd_rates[2395] - rev_rates[1253]);

  //rxn 2396
  //sp 114
  sp_rates[113] += 2.0 * (fwd_rates[2396] - rev_rates[1254]);
  //sp 115
  sp_rates[114] -= (fwd_rates[2396] - rev_rates[1254]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2396] - rev_rates[1254]);

  //rxn 2397
  //sp 115
  sp_rates[114] -= (fwd_rates[2397] - rev_rates[1255]);
  //sp 5
  sp_rates[4] += (fwd_rates[2397] - rev_rates[1255]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2397] - rev_rates[1255]);
  //sp 1
  (*dy_N) += (fwd_rates[2397] - rev_rates[1255]);

  //rxn 2398
  //sp 10
  sp_rates[9] += (fwd_rates[2398] - rev_rates[1256]);
  //sp 115
  sp_rates[114] -= (fwd_rates[2398] - rev_rates[1256]);
  //sp 1
  (*dy_N) += (fwd_rates[2398] - rev_rates[1256]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2398] - rev_rates[1256]);

  //rxn 2399
  //sp 114
  sp_rates[113] -= (fwd_rates[2399] - rev_rates[1257]);
  //sp 115
  sp_rates[114] -= (fwd_rates[2399] - rev_rates[1257]);
  //sp 116
  sp_rates[115] += (fwd_rates[2399] - rev_rates[1257]);
  //sp 1
  (*dy_N) += (fwd_rates[2399] - rev_rates[1257]);

  //rxn 2400
  //sp 114
  sp_rates[113] -= (fwd_rates[2400] - rev_rates[1258]);
  //sp 116
  sp_rates[115] += (fwd_rates[2400] - rev_rates[1258]);
  //sp 21
  sp_rates[20] -= (fwd_rates[2400] - rev_rates[1258]);
  //sp 23
  sp_rates[22] += (fwd_rates[2400] - rev_rates[1258]);

  //rxn 2401
  //sp 33
  sp_rates[32] -= (fwd_rates[2401] - rev_rates[1259]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2401] - rev_rates[1259]);
  //sp 116
  sp_rates[115] += (fwd_rates[2401] - rev_rates[1259]);
  //sp 39
  sp_rates[38] += (fwd_rates[2401] - rev_rates[1259]);

  //rxn 2402
  //sp 8
  sp_rates[7] += (fwd_rates[2402] - rev_rates[1260]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2402] - rev_rates[1260]);
  //sp 116
  sp_rates[115] += (fwd_rates[2402] - rev_rates[1260]);
  //sp 55
  sp_rates[54] -= (fwd_rates[2402] - rev_rates[1260]);
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[2402] - rev_rates[1260]);

  //rxn 2403
  //sp 114
  sp_rates[113] += (fwd_rates[2403] - rev_rates[1261]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2403] - rev_rates[1261]);
  //sp 15
  sp_rates[14] -= (fwd_rates[2403] - rev_rates[1261]);
  //sp 23
  sp_rates[22] += (fwd_rates[2403] - rev_rates[1261]);

  //rxn 2404
  //sp 39
  sp_rates[38] += (fwd_rates[2404] - rev_rates[1262]);
  //sp 114
  sp_rates[113] += (fwd_rates[2404] - rev_rates[1262]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2404] - rev_rates[1262]);
  //sp 31
  sp_rates[30] -= (fwd_rates[2404] - rev_rates[1262]);

  //rxn 2405
  //sp 35
  sp_rates[34] += fwd_rates[2405];
  //sp 114
  sp_rates[113] += fwd_rates[2405];
  //sp 116
  sp_rates[115] -= fwd_rates[2405];
  //sp 25
  sp_rates[24] += fwd_rates[2405];
  //sp 63
  sp_rates[62] -= fwd_rates[2405];

  //rxn 2406
  //sp 114
  sp_rates[113] += (fwd_rates[2406] - rev_rates[1263]);
  //sp 35
  sp_rates[34] -= (fwd_rates[2406] - rev_rates[1263]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2406] - rev_rates[1263]);
  //sp 47
  sp_rates[46] += (fwd_rates[2406] - rev_rates[1263]);

  //rxn 2407
  //sp 5
  sp_rates[4] += fwd_rates[2407];
  //sp 114
  sp_rates[113] += fwd_rates[2407];
  //sp 116
  sp_rates[115] -= fwd_rates[2407];
  //sp 21
  sp_rates[20] -= fwd_rates[2407];
  //sp 23
  sp_rates[22] += fwd_rates[2407];

  //rxn 2408
  //sp 33
  sp_rates[32] -= fwd_rates[2408];
  //sp 5
  sp_rates[4] += fwd_rates[2408];
  //sp 39
  sp_rates[38] += fwd_rates[2408];
  //sp 114
  sp_rates[113] += fwd_rates[2408];
  //sp 116
  sp_rates[115] -= fwd_rates[2408];

  //rxn 2409
  //sp 114
  sp_rates[113] -= (fwd_rates[2409] - rev_rates[1264]);
  //sp 117
  sp_rates[116] += (fwd_rates[2409] - rev_rates[1264]);
  //sp 30
  sp_rates[29] -= (fwd_rates[2409] - rev_rates[1264]);
  //sp 31
  sp_rates[30] += (fwd_rates[2409] - rev_rates[1264]);

  //rxn 2410
  //sp 4
  sp_rates[3] -= (fwd_rates[2410] - rev_rates[1265]);
  //sp 118
  sp_rates[117] += (fwd_rates[2410] - rev_rates[1265]);
  //sp 143
  sp_rates[142] = -(fwd_rates[2410] - rev_rates[1265]);
  //sp 15
  sp_rates[14] += (fwd_rates[2410] - rev_rates[1265]);

  //rxn 2411
  //sp 116
  sp_rates[115] += (fwd_rates[2411] - rev_rates[1266]) * pres_mod[60];
  //sp 143
  sp_rates[142] -= (fwd_rates[2411] - rev_rates[1266]) * pres_mod[60];
  //sp 15
  sp_rates[14] += (fwd_rates[2411] - rev_rates[1266]) * pres_mod[60];

  //rxn 2412
  //sp 118
  sp_rates[117] -= (fwd_rates[2412] - rev_rates[1267]);
  //sp 116
  sp_rates[115] += (fwd_rates[2412] - rev_rates[1267]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2412] - rev_rates[1267]);
  //sp 8
  sp_rates[7] += (fwd_rates[2412] - rev_rates[1267]);

  //rxn 2413
  //sp 116
  sp_rates[115] += (fwd_rates[2413] - rev_rates[1268]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2413] - rev_rates[1268]);
  //sp 7
  sp_rates[6] += (fwd_rates[2413] - rev_rates[1268]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2413] - rev_rates[1268]);

  //rxn 2414
  //sp 118
  sp_rates[117] -= (fwd_rates[2414] - rev_rates[1269]);
  //sp 119
  sp_rates[118] += (fwd_rates[2414] - rev_rates[1269]);

  //rxn 2415
  //sp 116
  sp_rates[115] += (fwd_rates[2415] - rev_rates[1270]);
  //sp 133
  sp_rates[132] -= (fwd_rates[2415] - rev_rates[1270]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2415] - rev_rates[1270]);
  //sp 127
  sp_rates[126] += (fwd_rates[2415] - rev_rates[1270]);

  //rxn 2416
  //sp 114
  sp_rates[113] += (fwd_rates[2416] - rev_rates[1271]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2416] - rev_rates[1271]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2416] - rev_rates[1271]);
  //sp 7
  sp_rates[6] += (fwd_rates[2416] - rev_rates[1271]);

  //rxn 2417
  //sp 4
  sp_rates[3] -= (fwd_rates[2417] - rev_rates[1272]);
  //sp 117
  sp_rates[116] += (fwd_rates[2417] - rev_rates[1272]);
  //sp 118
  sp_rates[117] -= (fwd_rates[2417] - rev_rates[1272]);
  //sp 8
  sp_rates[7] += (fwd_rates[2417] - rev_rates[1272]);

  //rxn 2418
  //sp 114
  sp_rates[113] -= (fwd_rates[2418] - rev_rates[1273]);
  //sp 118
  sp_rates[117] += (fwd_rates[2418] - rev_rates[1273]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2418] - rev_rates[1273]);

  //rxn 2419
  //sp 116
  sp_rates[115] -= (fwd_rates[2419] - rev_rates[1274]) * pres_mod[61];
  //sp 120
  sp_rates[119] = (fwd_rates[2419] - rev_rates[1274]) * pres_mod[61];
  //sp 8
  sp_rates[7] -= (fwd_rates[2419] - rev_rates[1274]) * pres_mod[61];

  //rxn 2420
  //sp 130
  sp_rates[129] += (fwd_rates[2420] - rev_rates[1275]);
  //sp 120
  sp_rates[119] -= (fwd_rates[2420] - rev_rates[1275]);
  //sp 7
  sp_rates[6] += (fwd_rates[2420] - rev_rates[1275]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2420] - rev_rates[1275]);

  //rxn 2421
  //sp 114
  sp_rates[113] -= (fwd_rates[2421] - rev_rates[1276]) * pres_mod[62];
  //sp 23
  sp_rates[22] -= (fwd_rates[2421] - rev_rates[1276]) * pres_mod[62];
  //sp 144
  sp_rates[143] = (fwd_rates[2421] - rev_rates[1276]) * pres_mod[62];

  //rxn 2422
  //sp 25
  sp_rates[24] += (fwd_rates[2422] - rev_rates[1277]);
  //sp 114
  sp_rates[113] -= (fwd_rates[2422] - rev_rates[1277]);
  //sp 117
  sp_rates[116] += (fwd_rates[2422] - rev_rates[1277]);
  //sp 23
  sp_rates[22] -= (fwd_rates[2422] - rev_rates[1277]);

  //rxn 2423
  //sp 114
  sp_rates[113] -= (fwd_rates[2423] - rev_rates[1278]) * pres_mod[63];
  //sp 142
  sp_rates[141] += (fwd_rates[2423] - rev_rates[1278]) * pres_mod[63];
  //sp 15
  sp_rates[14] -= (fwd_rates[2423] - rev_rates[1278]) * pres_mod[63];

  //rxn 2424
  //sp 145
  sp_rates[144] = (fwd_rates[2424] - rev_rates[1279]) * pres_mod[64];
  //sp 116
  sp_rates[115] -= (fwd_rates[2424] - rev_rates[1279]) * pres_mod[64];
  //sp 23
  sp_rates[22] -= (fwd_rates[2424] - rev_rates[1279]) * pres_mod[64];

  //rxn 2425
  //sp 145
  sp_rates[144] -= (fwd_rates[2425] - rev_rates[1280]);
  //sp 119
  sp_rates[118] += (fwd_rates[2425] - rev_rates[1280]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2425] - rev_rates[1280]);
  //sp 23
  sp_rates[22] += (fwd_rates[2425] - rev_rates[1280]);

  //rxn 2426
  //sp 145
  sp_rates[144] -= (fwd_rates[2426] - rev_rates[1281]);
  //sp 130
  sp_rates[129] += (fwd_rates[2426] - rev_rates[1281]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2426] - rev_rates[1281]);
  //sp 23
  sp_rates[22] += (fwd_rates[2426] - rev_rates[1281]);

  //rxn 2427
  //sp 145
  sp_rates[144] -= (fwd_rates[2427] - rev_rates[1282]);
  //sp 120
  sp_rates[119] += (fwd_rates[2427] - rev_rates[1282]);
  //sp 23
  sp_rates[22] += (fwd_rates[2427] - rev_rates[1282]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2427] - rev_rates[1282]);

  //rxn 2428
  //sp 25
  sp_rates[24] += (fwd_rates[2428] - rev_rates[1283]);
  //sp 119
  sp_rates[118] += (fwd_rates[2428] - rev_rates[1283]);
  //sp 116
  sp_rates[115] -= (fwd_rates[2428] - rev_rates[1283]);
  //sp 23
  sp_rates[22] -= (fwd_rates[2428] - rev_rates[1283]);

  //rxn 2429
  //sp 19
  sp_rates[18] -= (fwd_rates[2429] - rev_rates[1284]);
  //sp 4
  sp_rates[3] += (fwd_rates[2429] - rev_rates[1284]);
  //sp 148
  sp_rates[147] = (fwd_rates[2429] - rev_rates[1284]);
  //sp 1
  (*dy_N) -= (fwd_rates[2429] - rev_rates[1284]);

  //rxn 2430
  //sp 147
  sp_rates[146] -= (fwd_rates[2430] - rev_rates[1285]);
  //sp 11
  sp_rates[10] += (fwd_rates[2430] - rev_rates[1285]);
  //sp 149
  sp_rates[148] -= (fwd_rates[2430] - rev_rates[1285]);
  //sp 148
  sp_rates[147] += (fwd_rates[2430] - rev_rates[1285]);

  //rxn 2431
  //sp 19
  sp_rates[18] -= (fwd_rates[2431] - rev_rates[1286]);
  //sp 150
  sp_rates[149] = (fwd_rates[2431] - rev_rates[1286]);
  //sp 1
  (*dy_N) -= (fwd_rates[2431] - rev_rates[1286]);

  //rxn 2432
  //sp 131
  sp_rates[130] += (fwd_rates[2432] - rev_rates[1287]);
  //sp 149
  sp_rates[148] += (fwd_rates[2432] - rev_rates[1287]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2432] - rev_rates[1287]);
  //sp 150
  sp_rates[149] -= (fwd_rates[2432] - rev_rates[1287]);

  //rxn 2433
  //sp 147
  sp_rates[146] += (fwd_rates[2433] - rev_rates[1288]);
  //sp 117
  sp_rates[116] += (fwd_rates[2433] - rev_rates[1288]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2433] - rev_rates[1288]);
  //sp 150
  sp_rates[149] -= (fwd_rates[2433] - rev_rates[1288]);

  //rxn 2434
  //sp 148
  sp_rates[147] += (fwd_rates[2434] - rev_rates[1289]);
  //sp 150
  sp_rates[149] -= (fwd_rates[2434] - rev_rates[1289]);
  //sp 7
  sp_rates[6] += (fwd_rates[2434] - rev_rates[1289]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2434] - rev_rates[1289]);

  //rxn 2435
  //sp 10
  sp_rates[9] += (fwd_rates[2435] - rev_rates[1290]);
  //sp 148
  sp_rates[147] += (fwd_rates[2435] - rev_rates[1290]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2435] - rev_rates[1290]);
  //sp 150
  sp_rates[149] -= (fwd_rates[2435] - rev_rates[1290]);

  //rxn 2436
  //sp 129
  sp_rates[128] += (fwd_rates[2436] - rev_rates[1291]);
  //sp 147
  sp_rates[146] += (fwd_rates[2436] - rev_rates[1291]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2436] - rev_rates[1291]);

  //rxn 2437
  //sp 18
  sp_rates[17] += (fwd_rates[2437] - rev_rates[1292]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2437] - rev_rates[1292]);
  //sp 1
  (*dy_N) += (fwd_rates[2437] - rev_rates[1292]);

  //rxn 2438
  //sp 129
  sp_rates[128] += (fwd_rates[2438] - rev_rates[1293]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2438] - rev_rates[1293]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2438] - rev_rates[1293]);
  //sp 136
  sp_rates[135] += (fwd_rates[2438] - rev_rates[1293]);

  //rxn 2439
  //sp 4
  sp_rates[3] -= (fwd_rates[2439] - rev_rates[1294]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2439] - rev_rates[1294]);
  //sp 150
  sp_rates[149] += (fwd_rates[2439] - rev_rates[1294]);

  //rxn 2440
  //sp 114
  sp_rates[113] += (fwd_rates[2440] - rev_rates[1295]);
  //sp 147
  sp_rates[146] += (fwd_rates[2440] - rev_rates[1295]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2440] - rev_rates[1295]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2440] - rev_rates[1295]);

  //rxn 2441
  //sp 11
  sp_rates[10] += (fwd_rates[2441] - rev_rates[1296]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2441] - rev_rates[1296]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2441] - rev_rates[1296]);
  //sp 1
  (*dy_N) += (fwd_rates[2441] - rev_rates[1296]);

  //rxn 2442
  //sp 129
  sp_rates[128] += (fwd_rates[2442] - rev_rates[1297]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2442] - rev_rates[1297]);
  //sp 149
  sp_rates[148] += (fwd_rates[2442] - rev_rates[1297]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2442] - rev_rates[1297]);

  //rxn 2443
  //sp 136
  sp_rates[135] += (fwd_rates[2443] - rev_rates[1298]);
  //sp 114
  sp_rates[113] += (fwd_rates[2443] - rev_rates[1298]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2443] - rev_rates[1298]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2443] - rev_rates[1298]);

  //rxn 2444
  //sp 114
  sp_rates[113] += (fwd_rates[2444] - rev_rates[1299]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2444] - rev_rates[1299]);
  //sp 5
  sp_rates[4] -= (fwd_rates[2444] - rev_rates[1299]);
  //sp 149
  sp_rates[148] += (fwd_rates[2444] - rev_rates[1299]);

  //rxn 2445
  //sp 19
  sp_rates[18] -= (fwd_rates[2445] - rev_rates[1300]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2445] - rev_rates[1300]);
  //sp 147
  sp_rates[146] += (fwd_rates[2445] - rev_rates[1300]);
  //sp 136
  sp_rates[135] += (fwd_rates[2445] - rev_rates[1300]);

  //rxn 2446
  //sp 147
  sp_rates[146] += (fwd_rates[2446] - rev_rates[1301]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2446] - rev_rates[1301]);
  //sp 151
  sp_rates[150] += (fwd_rates[2446] - rev_rates[1301]);
  //sp 16
  sp_rates[15] -= (fwd_rates[2446] - rev_rates[1301]);

  //rxn 2447
  //sp 4
  sp_rates[3] += (fwd_rates[2447] - rev_rates[1302]);
  //sp 3
  sp_rates[2] -= (fwd_rates[2447] - rev_rates[1302]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2447] - rev_rates[1302]);
  //sp 150
  sp_rates[149] += (fwd_rates[2447] - rev_rates[1302]);

  //rxn 2448
  //sp 18
  sp_rates[17] += (fwd_rates[2448] - rev_rates[1303]) * pres_mod[65];
  //sp 148
  sp_rates[147] -= (fwd_rates[2448] - rev_rates[1303]) * pres_mod[65];
  //sp 1
  (*dy_N) += (fwd_rates[2448] - rev_rates[1303]) * pres_mod[65];

  //rxn 2449
  //sp 147
  sp_rates[146] += 2.0 * (fwd_rates[2449] - rev_rates[1304]);
  //sp 148
  sp_rates[147] -= 2.0 * (fwd_rates[2449] - rev_rates[1304]);
  //sp 1
  (*dy_N) += (fwd_rates[2449] - rev_rates[1304]);

  //rxn 2450
  //sp 18
  sp_rates[17] -= (fwd_rates[2450] - rev_rates[1305]);
  //sp 147
  sp_rates[146] += 2.0 * (fwd_rates[2450] - rev_rates[1305]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2450] - rev_rates[1305]);

  //rxn 2451
  //sp 129
  sp_rates[128] -= (fwd_rates[2451] - rev_rates[1306]);
  //sp 147
  sp_rates[146] += (fwd_rates[2451] - rev_rates[1306]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2451] - rev_rates[1306]);
  //sp 1
  (*dy_N) += (fwd_rates[2451] - rev_rates[1306]);

  //rxn 2452
  //sp 129
  sp_rates[128] += (fwd_rates[2452] - rev_rates[1307]);
  //sp 147
  sp_rates[146] -= (fwd_rates[2452] - rev_rates[1307]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2452] - rev_rates[1307]);
  //sp 153
  sp_rates[152] += (fwd_rates[2452] - rev_rates[1307]);

  //rxn 2453
  //sp 114
  sp_rates[113] -= (fwd_rates[2453] - rev_rates[1308]);
  //sp 115
  sp_rates[114] += (fwd_rates[2453] - rev_rates[1308]);
  //sp 148
  sp_rates[147] -= (fwd_rates[2453] - rev_rates[1308]);
  //sp 147
  sp_rates[146] += (fwd_rates[2453] - rev_rates[1308]);

  //rxn 2454
  //sp 153
  sp_rates[152] -= (fwd_rates[2454] - rev_rates[1309]);
  //sp 147
  sp_rates[146] += (fwd_rates[2454] - rev_rates[1309]);
  //sp 149
  sp_rates[148] += (fwd_rates[2454] - rev_rates[1309]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2454] - rev_rates[1309]);

  //rxn 2455
  //sp 153
  sp_rates[152] -= (fwd_rates[2455] - rev_rates[1310]);
  //sp 147
  sp_rates[146] += (fwd_rates[2455] - rev_rates[1310]);
  //sp 140
  sp_rates[139] += (fwd_rates[2455] - rev_rates[1310]);
  //sp 8
  sp_rates[7] -= (fwd_rates[2455] - rev_rates[1310]);

  //rxn 2456
  //sp 153
  sp_rates[152] -= (fwd_rates[2456] - rev_rates[1311]) * pres_mod[66];
  //sp 147
  sp_rates[146] += 2.0 * (fwd_rates[2456] - rev_rates[1311]) * pres_mod[66];

  //rxn 2457
  //sp 153
  sp_rates[152] -= (fwd_rates[2457] - rev_rates[1312]);
  //sp 147
  sp_rates[146] += (fwd_rates[2457] - rev_rates[1312]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2457] - rev_rates[1312]);
  //sp 136
  sp_rates[135] += (fwd_rates[2457] - rev_rates[1312]);

  //rxn 2458
  //sp 147
  sp_rates[146] -= (fwd_rates[2458] - rev_rates[1313]);
  //sp 36
  sp_rates[35] -= (fwd_rates[2458] - rev_rates[1313]);
  //sp 37
  sp_rates[36] += (fwd_rates[2458] - rev_rates[1313]);
  //sp 136
  sp_rates[135] += (fwd_rates[2458] - rev_rates[1313]);

  //sp 0
  sp_rates[0] = 0.0;
  //sp 2
  sp_rates[1] = 0.0;
  //sp 52
  sp_rates[51] = 0.0;
} // end eval_spec_rates

