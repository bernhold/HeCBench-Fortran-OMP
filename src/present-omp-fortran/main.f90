program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(int32), save :: sbox(0:15) = [ &
    192_int32, 80_int32, 96_int32, 176_int32, 144_int32, 0_int32, 160_int32, 208_int32, 48_int32, 224_int32, 240_int32, 128_int32, 64_int32, 112_int32, 16_int32, 32_int32 &
  ]
  integer(int32), save :: sbox_pmt_3(0:255) = [ &
    240_int32, 177_int32, 180_int32, 229_int32, 225_int32, 160_int32, 228_int32, 241_int32, 165_int32, 244_int32, 245_int32, 224_int32, 176_int32, 181_int32, 161_int32, 164_int32, &
    114_int32, 51_int32, 54_int32, 103_int32, 99_int32, 34_int32, 102_int32, 115_int32, 39_int32, 118_int32, 119_int32, 98_int32, 50_int32, 55_int32, 35_int32, 38_int32, &
    120_int32, 57_int32, 60_int32, 109_int32, 105_int32, 40_int32, 108_int32, 121_int32, 45_int32, 124_int32, 125_int32, 104_int32, 56_int32, 61_int32, 41_int32, 44_int32, &
    218_int32, 155_int32, 158_int32, 207_int32, 203_int32, 138_int32, 206_int32, 219_int32, 143_int32, 222_int32, 223_int32, 202_int32, 154_int32, 159_int32, 139_int32, 142_int32, &
    210_int32, 147_int32, 150_int32, 199_int32, 195_int32, 130_int32, 198_int32, 211_int32, 135_int32, 214_int32, 215_int32, 194_int32, 146_int32, 151_int32, 131_int32, 134_int32, &
    80_int32, 17_int32, 20_int32, 69_int32, 65_int32, 0_int32, 68_int32, 81_int32, 5_int32, 84_int32, 85_int32, 64_int32, 16_int32, 21_int32, 1_int32, 4_int32, &
    216_int32, 153_int32, 156_int32, 205_int32, 201_int32, 136_int32, 204_int32, 217_int32, 141_int32, 220_int32, 221_int32, 200_int32, 152_int32, 157_int32, 137_int32, 140_int32, &
    242_int32, 179_int32, 182_int32, 231_int32, 227_int32, 162_int32, 230_int32, 243_int32, 167_int32, 246_int32, 247_int32, 226_int32, 178_int32, 183_int32, 163_int32, 166_int32, &
    90_int32, 27_int32, 30_int32, 79_int32, 75_int32, 10_int32, 78_int32, 91_int32, 15_int32, 94_int32, 95_int32, 74_int32, 26_int32, 31_int32, 11_int32, 14_int32, &
    248_int32, 185_int32, 188_int32, 237_int32, 233_int32, 168_int32, 236_int32, 249_int32, 173_int32, 252_int32, 253_int32, 232_int32, 184_int32, 189_int32, 169_int32, 172_int32, &
    250_int32, 187_int32, 190_int32, 239_int32, 235_int32, 170_int32, 238_int32, 251_int32, 175_int32, 254_int32, 255_int32, 234_int32, 186_int32, 191_int32, 171_int32, 174_int32, &
    208_int32, 145_int32, 148_int32, 197_int32, 193_int32, 128_int32, 196_int32, 209_int32, 133_int32, 212_int32, 213_int32, 192_int32, 144_int32, 149_int32, 129_int32, 132_int32, &
    112_int32, 49_int32, 52_int32, 101_int32, 97_int32, 32_int32, 100_int32, 113_int32, 37_int32, 116_int32, 117_int32, 96_int32, 48_int32, 53_int32, 33_int32, 36_int32, &
    122_int32, 59_int32, 62_int32, 111_int32, 107_int32, 42_int32, 110_int32, 123_int32, 47_int32, 126_int32, 127_int32, 106_int32, 58_int32, 63_int32, 43_int32, 46_int32, &
    82_int32, 19_int32, 22_int32, 71_int32, 67_int32, 2_int32, 70_int32, 83_int32, 7_int32, 86_int32, 87_int32, 66_int32, 18_int32, 23_int32, 3_int32, 6_int32, &
    88_int32, 25_int32, 28_int32, 77_int32, 73_int32, 8_int32, 76_int32, 89_int32, 13_int32, 92_int32, 93_int32, 72_int32, 24_int32, 29_int32, 9_int32, 12_int32 &
  ]
  integer(int32), save :: sbox_pmt_2(0:255) = [ &
    60_int32, 108_int32, 45_int32, 121_int32, 120_int32, 40_int32, 57_int32, 124_int32, 105_int32, 61_int32, 125_int32, 56_int32, 44_int32, 109_int32, 104_int32, 41_int32, &
    156_int32, 204_int32, 141_int32, 217_int32, 216_int32, 136_int32, 153_int32, 220_int32, 201_int32, 157_int32, 221_int32, 152_int32, 140_int32, 205_int32, 200_int32, 137_int32, &
    30_int32, 78_int32, 15_int32, 91_int32, 90_int32, 10_int32, 27_int32, 94_int32, 75_int32, 31_int32, 95_int32, 26_int32, 14_int32, 79_int32, 74_int32, 11_int32, &
    182_int32, 230_int32, 167_int32, 243_int32, 242_int32, 162_int32, 179_int32, 246_int32, 227_int32, 183_int32, 247_int32, 178_int32, 166_int32, 231_int32, 226_int32, 163_int32, &
    180_int32, 228_int32, 165_int32, 241_int32, 240_int32, 160_int32, 177_int32, 244_int32, 225_int32, 181_int32, 245_int32, 176_int32, 164_int32, 229_int32, 224_int32, 161_int32, &
    20_int32, 68_int32, 5_int32, 81_int32, 80_int32, 0_int32, 17_int32, 84_int32, 65_int32, 21_int32, 85_int32, 16_int32, 4_int32, 69_int32, 64_int32, 1_int32, &
    54_int32, 102_int32, 39_int32, 115_int32, 114_int32, 34_int32, 51_int32, 118_int32, 99_int32, 55_int32, 119_int32, 50_int32, 38_int32, 103_int32, 98_int32, 35_int32, &
    188_int32, 236_int32, 173_int32, 249_int32, 248_int32, 168_int32, 185_int32, 252_int32, 233_int32, 189_int32, 253_int32, 184_int32, 172_int32, 237_int32, 232_int32, 169_int32, &
    150_int32, 198_int32, 135_int32, 211_int32, 210_int32, 130_int32, 147_int32, 214_int32, 195_int32, 151_int32, 215_int32, 146_int32, 134_int32, 199_int32, 194_int32, 131_int32, &
    62_int32, 110_int32, 47_int32, 123_int32, 122_int32, 42_int32, 59_int32, 126_int32, 107_int32, 63_int32, 127_int32, 58_int32, 46_int32, 111_int32, 106_int32, 43_int32, &
    190_int32, 238_int32, 175_int32, 251_int32, 250_int32, 170_int32, 187_int32, 254_int32, 235_int32, 191_int32, 255_int32, 186_int32, 174_int32, 239_int32, 234_int32, 171_int32, &
    52_int32, 100_int32, 37_int32, 113_int32, 112_int32, 32_int32, 49_int32, 116_int32, 97_int32, 53_int32, 117_int32, 48_int32, 36_int32, 101_int32, 96_int32, 33_int32, &
    28_int32, 76_int32, 13_int32, 89_int32, 88_int32, 8_int32, 25_int32, 92_int32, 73_int32, 29_int32, 93_int32, 24_int32, 12_int32, 77_int32, 72_int32, 9_int32, &
    158_int32, 206_int32, 143_int32, 219_int32, 218_int32, 138_int32, 155_int32, 222_int32, 203_int32, 159_int32, 223_int32, 154_int32, 142_int32, 207_int32, 202_int32, 139_int32, &
    148_int32, 196_int32, 133_int32, 209_int32, 208_int32, 128_int32, 145_int32, 212_int32, 193_int32, 149_int32, 213_int32, 144_int32, 132_int32, 197_int32, 192_int32, 129_int32, &
    22_int32, 70_int32, 7_int32, 83_int32, 82_int32, 2_int32, 19_int32, 86_int32, 67_int32, 23_int32, 87_int32, 18_int32, 6_int32, 71_int32, 66_int32, 3_int32 &
  ]
  integer(int32), save :: sbox_pmt_1(0:255) = [ &
    15_int32, 27_int32, 75_int32, 94_int32, 30_int32, 10_int32, 78_int32, 31_int32, 90_int32, 79_int32, 95_int32, 14_int32, 11_int32, 91_int32, 26_int32, 74_int32, &
    39_int32, 51_int32, 99_int32, 118_int32, 54_int32, 34_int32, 102_int32, 55_int32, 114_int32, 103_int32, 119_int32, 38_int32, 35_int32, 115_int32, 50_int32, 98_int32, &
    135_int32, 147_int32, 195_int32, 214_int32, 150_int32, 130_int32, 198_int32, 151_int32, 210_int32, 199_int32, 215_int32, 134_int32, 131_int32, 211_int32, 146_int32, 194_int32, &
    173_int32, 185_int32, 233_int32, 252_int32, 188_int32, 168_int32, 236_int32, 189_int32, 248_int32, 237_int32, 253_int32, 172_int32, 169_int32, 249_int32, 184_int32, 232_int32, &
    45_int32, 57_int32, 105_int32, 124_int32, 60_int32, 40_int32, 108_int32, 61_int32, 120_int32, 109_int32, 125_int32, 44_int32, 41_int32, 121_int32, 56_int32, 104_int32, &
    5_int32, 17_int32, 65_int32, 84_int32, 20_int32, 0_int32, 68_int32, 21_int32, 80_int32, 69_int32, 85_int32, 4_int32, 1_int32, 81_int32, 16_int32, 64_int32, &
    141_int32, 153_int32, 201_int32, 220_int32, 156_int32, 136_int32, 204_int32, 157_int32, 216_int32, 205_int32, 221_int32, 140_int32, 137_int32, 217_int32, 152_int32, 200_int32, &
    47_int32, 59_int32, 107_int32, 126_int32, 62_int32, 42_int32, 110_int32, 63_int32, 122_int32, 111_int32, 127_int32, 46_int32, 43_int32, 123_int32, 58_int32, 106_int32, &
    165_int32, 177_int32, 225_int32, 244_int32, 180_int32, 160_int32, 228_int32, 181_int32, 240_int32, 229_int32, 245_int32, 164_int32, 161_int32, 241_int32, 176_int32, 224_int32, &
    143_int32, 155_int32, 203_int32, 222_int32, 158_int32, 138_int32, 206_int32, 159_int32, 218_int32, 207_int32, 223_int32, 142_int32, 139_int32, 219_int32, 154_int32, 202_int32, &
    175_int32, 187_int32, 235_int32, 254_int32, 190_int32, 170_int32, 238_int32, 191_int32, 250_int32, 239_int32, 255_int32, 174_int32, 171_int32, 251_int32, 186_int32, 234_int32, &
    13_int32, 25_int32, 73_int32, 92_int32, 28_int32, 8_int32, 76_int32, 29_int32, 88_int32, 77_int32, 93_int32, 12_int32, 9_int32, 89_int32, 24_int32, 72_int32, &
    7_int32, 19_int32, 67_int32, 86_int32, 22_int32, 2_int32, 70_int32, 23_int32, 82_int32, 71_int32, 87_int32, 6_int32, 3_int32, 83_int32, 18_int32, 66_int32, &
    167_int32, 179_int32, 227_int32, 246_int32, 182_int32, 162_int32, 230_int32, 183_int32, 242_int32, 231_int32, 247_int32, 166_int32, 163_int32, 243_int32, 178_int32, 226_int32, &
    37_int32, 49_int32, 97_int32, 116_int32, 52_int32, 32_int32, 100_int32, 53_int32, 112_int32, 101_int32, 117_int32, 36_int32, 33_int32, 113_int32, 48_int32, 96_int32, &
    133_int32, 145_int32, 193_int32, 212_int32, 148_int32, 128_int32, 196_int32, 149_int32, 208_int32, 197_int32, 213_int32, 132_int32, 129_int32, 209_int32, 144_int32, 192_int32 &
  ]
  integer(int32), save :: sbox_pmt_0(0:255) = [ &
    195_int32, 198_int32, 210_int32, 151_int32, 135_int32, 130_int32, 147_int32, 199_int32, 150_int32, 211_int32, 215_int32, 131_int32, 194_int32, 214_int32, 134_int32, 146_int32, &
    201_int32, 204_int32, 216_int32, 157_int32, 141_int32, 136_int32, 153_int32, 205_int32, 156_int32, 217_int32, 221_int32, 137_int32, 200_int32, 220_int32, 140_int32, 152_int32, &
    225_int32, 228_int32, 240_int32, 181_int32, 165_int32, 160_int32, 177_int32, 229_int32, 180_int32, 241_int32, 245_int32, 161_int32, 224_int32, 244_int32, 164_int32, 176_int32, &
    107_int32, 110_int32, 122_int32, 63_int32, 47_int32, 42_int32, 59_int32, 111_int32, 62_int32, 123_int32, 127_int32, 43_int32, 106_int32, 126_int32, 46_int32, 58_int32, &
    75_int32, 78_int32, 90_int32, 31_int32, 15_int32, 10_int32, 27_int32, 79_int32, 30_int32, 91_int32, 95_int32, 11_int32, 74_int32, 94_int32, 14_int32, 26_int32, &
    65_int32, 68_int32, 80_int32, 21_int32, 5_int32, 0_int32, 17_int32, 69_int32, 20_int32, 81_int32, 85_int32, 1_int32, 64_int32, 84_int32, 4_int32, 16_int32, &
    99_int32, 102_int32, 114_int32, 55_int32, 39_int32, 34_int32, 51_int32, 103_int32, 54_int32, 115_int32, 119_int32, 35_int32, 98_int32, 118_int32, 38_int32, 50_int32, &
    203_int32, 206_int32, 218_int32, 159_int32, 143_int32, 138_int32, 155_int32, 207_int32, 158_int32, 219_int32, 223_int32, 139_int32, 202_int32, 222_int32, 142_int32, 154_int32, &
    105_int32, 108_int32, 120_int32, 61_int32, 45_int32, 40_int32, 57_int32, 109_int32, 60_int32, 121_int32, 125_int32, 41_int32, 104_int32, 124_int32, 44_int32, 56_int32, &
    227_int32, 230_int32, 242_int32, 183_int32, 167_int32, 162_int32, 179_int32, 231_int32, 182_int32, 243_int32, 247_int32, 163_int32, 226_int32, 246_int32, 166_int32, 178_int32, &
    235_int32, 238_int32, 250_int32, 191_int32, 175_int32, 170_int32, 187_int32, 239_int32, 190_int32, 251_int32, 255_int32, 171_int32, 234_int32, 254_int32, 174_int32, 186_int32, &
    67_int32, 70_int32, 82_int32, 23_int32, 7_int32, 2_int32, 19_int32, 71_int32, 22_int32, 83_int32, 87_int32, 3_int32, 66_int32, 86_int32, 6_int32, 18_int32, &
    193_int32, 196_int32, 208_int32, 149_int32, 133_int32, 128_int32, 145_int32, 197_int32, 148_int32, 209_int32, 213_int32, 129_int32, 192_int32, 212_int32, 132_int32, 144_int32, &
    233_int32, 236_int32, 248_int32, 189_int32, 173_int32, 168_int32, 185_int32, 237_int32, 188_int32, 249_int32, 253_int32, 169_int32, 232_int32, 252_int32, 172_int32, 184_int32, &
    73_int32, 76_int32, 88_int32, 29_int32, 13_int32, 8_int32, 25_int32, 77_int32, 28_int32, 89_int32, 93_int32, 9_int32, 72_int32, 92_int32, 12_int32, 24_int32, &
    97_int32, 100_int32, 112_int32, 53_int32, 37_int32, 32_int32, 49_int32, 101_int32, 52_int32, 113_int32, 117_int32, 33_int32, 96_int32, 116_int32, 36_int32, 48_int32 &
  ]
  !$omp declare target(sbox, sbox_pmt_3, sbox_pmt_2, sbox_pmt_1, sbox_pmt_0)


  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer :: num, repeat, argc
  character(len=256) :: arg
  integer(int32), allocatable :: h_plain(:), h_key(:), h_cipher(:), ciphers(:)
  integer(int64) :: h_checksum, d_checksum
  integer(int32) :: plain(0:7), key(0:9)
  integer :: i, k, n
  integer(int32), parameter :: rounds = 31_int32
  integer(int32), parameter :: seed = 8_int32
  real(real64) :: start_time, end_time, elapsed

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, arg)
    write(*,'(A,A,A)') 'Usage: ', trim(arg), ' <number of plain texts> <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) num
  call get_command_argument(2, arg)
  read(arg, *) repeat
  if (num <= 0 .or. repeat <= 0) stop 1

  allocate(h_plain(0:8 * num - 1), h_key(0:10 * num - 1), h_cipher(0:8 * num - 1))
  allocate(ciphers(0:8 * num - 1))

  call c_srand(int(seed, c_int))
  plain = [int(iachar('P'), int32), int(iachar('R'), int32), int(iachar('E'), int32), &
    int(iachar('S'), int32), int(iachar('E'), int32), int(iachar('N'), int32), &
    int(iachar('T'), int32), 0_int32]

  do i = 0, num - 1
    do k = 0, 9
      key(k) = mod(c_rand(), 256)
      h_key(i * 10 + k) = key(k)
    end do
    do k = 0, 7
      h_plain(i * 8 + k) = plain(k)
    end do
    call deterministic_shuffle(plain)
  end do

  h_checksum = 0_int64
  do n = 0, repeat
    do i = 0, num - 1
      call present_rounds(h_plain(i * 8:), h_key(i * 10:), rounds, h_cipher(i * 8:))
      do k = 0, 7
        h_checksum = h_checksum + int(h_cipher(i * 8 + k), int64)
      end do
    end do
  end do

  ciphers = 0_int32
  d_checksum = 0_int64
  elapsed = 0.0_real64

  !$omp target data map(to: h_plain(0:8 * num - 1), h_key(0:10 * num - 1)) &
  !$omp& map(alloc: ciphers(0:8 * num - 1))
  do n = 0, repeat
    start_time = omp_get_wtime()
    call present_kernel(num, h_plain, h_key, rounds, ciphers)
    end_time = omp_get_wtime()
    if (n > 0) elapsed = elapsed + (end_time - start_time)

    !$omp target update from(ciphers(0:8 * num - 1))
    do i = 0, 8 * num - 1
      d_checksum = d_checksum + int(ciphers(i), int64)
    end do
  end do
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', (elapsed * 1.0e6_real64) / real(repeat, real64), ' (us)'
  if (h_checksum /= d_checksum) then
    write(*,'(A)') 'FAIL'
  else
    write(*,'(A)') 'PASS'
  end if

  deallocate(ciphers, h_cipher, h_key, h_plain)

contains

  subroutine deterministic_shuffle(values)
    integer(int32), intent(inout) :: values(0:7)
    integer(int32), parameter :: order(0:7) = [4_int32, 5_int32, 1_int32, 2_int32, 6_int32, 7_int32, 3_int32, 0_int32]
    integer(int32) :: tmp(0:7)
    integer :: j

    tmp = values
    do j = 0, 7
      values(j) = tmp(order(j))
    end do
  end subroutine deterministic_shuffle

  subroutine present_kernel(num_items, plains, keys, rounds, ciphers)
    integer, intent(in) :: num_items
    integer(int32), intent(in) :: plains(0:), keys(0:)
    integer(int32), intent(in) :: rounds
    integer(int32), intent(inout) :: ciphers(0:)
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256)
    do idx = 0, num_items - 1
      call present_rounds(plains(idx * 8:), keys(idx * 10:), rounds, ciphers(idx * 8:))
    end do
  end subroutine present_kernel

  subroutine present_rounds(plain, key, rounds, cipher)
    !$omp declare target
    integer(int32), intent(in) :: plain(0:), key(0:)
    integer(int32), intent(in) :: rounds
    integer(int32), intent(inout) :: cipher(0:)
    integer(int32) :: round_counter
    integer(int32) :: state(0:7), round_key(0:9)

    state(0) = bxor(plain(0), key(0))
    state(1) = bxor(plain(1), key(1))
    state(2) = bxor(plain(2), key(2))
    state(3) = bxor(plain(3), key(3))
    state(4) = bxor(plain(4), key(4))
    state(5) = bxor(plain(5), key(5))
    state(6) = bxor(plain(6), key(6))
    state(7) = bxor(plain(7), key(7))

    round_key(9) = bor(shl8(key(6), 5), shr8(key(7), 3))
    round_key(8) = bor(shl8(key(5), 5), shr8(key(6), 3))
    round_key(7) = bor(shl8(key(4), 5), shr8(key(5), 3))
    round_key(6) = bor(shl8(key(3), 5), shr8(key(4), 3))
    round_key(5) = bor(shl8(key(2), 5), shr8(key(3), 3))
    round_key(4) = bor(shl8(key(1), 5), shr8(key(2), 3))
    round_key(3) = bor(shl8(key(0), 5), shr8(key(1), 3))
    round_key(2) = bor(shl8(key(9), 5), shr8(key(0), 3))
    round_key(1) = bor(shl8(key(8), 5), shr8(key(9), 3))
    round_key(0) = bor(shl8(key(7), 5), shr8(key(8), 3))
    round_key(0) = bor(band(round_key(0), z'0F'), sbox(shr8(round_key(0), 4)))
    round_key(7) = bxor(round_key(7), shr8(1_int32, 1))
    round_key(8) = bxor(round_key(8), shl8(1_int32, 7))

    call substitute_permute(state, cipher)

    do round_counter = 2_int32, rounds
      state(0) = bxor(cipher(0), round_key(0))
      state(1) = bxor(cipher(1), round_key(1))
      state(2) = bxor(cipher(2), round_key(2))
      state(3) = bxor(cipher(3), round_key(3))
      state(4) = bxor(cipher(4), round_key(4))
      state(5) = bxor(cipher(5), round_key(5))
      state(6) = bxor(cipher(6), round_key(6))
      state(7) = bxor(cipher(7), round_key(7))

      call substitute_permute(state, cipher)
      call update_round_key(round_key, state, round_counter)
    end do

    if (rounds == 31_int32) then
      cipher(0) = bxor(cipher(0), round_key(0))
      cipher(1) = bxor(cipher(1), round_key(1))
      cipher(2) = bxor(cipher(2), round_key(2))
      cipher(3) = bxor(cipher(3), round_key(3))
      cipher(4) = bxor(cipher(4), round_key(4))
      cipher(5) = bxor(cipher(5), round_key(5))
      cipher(6) = bxor(cipher(6), round_key(6))
      cipher(7) = bxor(cipher(7), round_key(7))
    end if
  end subroutine present_rounds

  subroutine substitute_permute(state, cipher)
    !$omp declare target
    integer(int32), intent(in) :: state(0:7)
    integer(int32), intent(inout) :: cipher(0:7)

    cipher(0) = bor(bor(bor(band(sbox_pmt_3(state(0)), z'C0'), band(sbox_pmt_2(state(1)), z'30')), &
      band(sbox_pmt_1(state(2)), z'0C')), band(sbox_pmt_0(state(3)), z'03'))
    cipher(1) = bor(bor(bor(band(sbox_pmt_3(state(4)), z'C0'), band(sbox_pmt_2(state(5)), z'30')), &
      band(sbox_pmt_1(state(6)), z'0C')), band(sbox_pmt_0(state(7)), z'03'))
    cipher(2) = bor(bor(bor(band(sbox_pmt_0(state(0)), z'C0'), band(sbox_pmt_3(state(1)), z'30')), &
      band(sbox_pmt_2(state(2)), z'0C')), band(sbox_pmt_1(state(3)), z'03'))
    cipher(3) = bor(bor(bor(band(sbox_pmt_0(state(4)), z'C0'), band(sbox_pmt_3(state(5)), z'30')), &
      band(sbox_pmt_2(state(6)), z'0C')), band(sbox_pmt_1(state(7)), z'03'))
    cipher(4) = bor(bor(bor(band(sbox_pmt_1(state(0)), z'C0'), band(sbox_pmt_0(state(1)), z'30')), &
      band(sbox_pmt_3(state(2)), z'0C')), band(sbox_pmt_2(state(3)), z'03'))
    cipher(5) = bor(bor(bor(band(sbox_pmt_1(state(4)), z'C0'), band(sbox_pmt_0(state(5)), z'30')), &
      band(sbox_pmt_3(state(6)), z'0C')), band(sbox_pmt_2(state(7)), z'03'))
    cipher(6) = bor(bor(bor(band(sbox_pmt_2(state(0)), z'C0'), band(sbox_pmt_1(state(1)), z'30')), &
      band(sbox_pmt_0(state(2)), z'0C')), band(sbox_pmt_3(state(3)), z'03'))
    cipher(7) = bor(bor(bor(band(sbox_pmt_2(state(4)), z'C0'), band(sbox_pmt_1(state(5)), z'30')), &
      band(sbox_pmt_0(state(6)), z'0C')), band(sbox_pmt_3(state(7)), z'03'))
  end subroutine substitute_permute

  subroutine update_round_key(round_key, scratch, round_counter)
    !$omp declare target
    integer(int32), intent(inout) :: round_key(0:9)
    integer(int32), intent(inout) :: scratch(0:7)
    integer(int32), intent(in) :: round_counter

    round_key(5) = bxor(round_key(5), shl8(round_counter, 2))
    scratch(2) = round_key(9)
    scratch(1) = round_key(8)
    scratch(0) = round_key(7)
    round_key(9) = bor(shl8(round_key(6), 5), shr8(round_key(7), 3))
    round_key(8) = bor(shl8(round_key(5), 5), shr8(round_key(6), 3))
    round_key(7) = bor(shl8(round_key(4), 5), shr8(round_key(5), 3))
    round_key(6) = bor(shl8(round_key(3), 5), shr8(round_key(4), 3))
    round_key(5) = bor(shl8(round_key(2), 5), shr8(round_key(3), 3))
    round_key(4) = bor(shl8(round_key(1), 5), shr8(round_key(2), 3))
    round_key(3) = bor(shl8(round_key(0), 5), shr8(round_key(1), 3))
    round_key(2) = bor(shl8(scratch(2), 5), shr8(round_key(0), 3))
    round_key(1) = bor(shl8(scratch(1), 5), shr8(scratch(2), 3))
    round_key(0) = bor(shl8(scratch(0), 5), shr8(scratch(1), 3))
    round_key(0) = bor(band(round_key(0), z'0F'), sbox(shr8(round_key(0), 4)))
  end subroutine update_round_key

  integer(int32) function band(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(lhs, rhs)
  end function band

  integer(int32) function bor(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(ior(lhs, rhs), 255_int32)
  end function bor

  integer(int32) function bxor(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(ieor(lhs, rhs), 255_int32)
  end function bxor

  integer(int32) function shl8(lhs, amount) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, amount
    value = iand(ishft(iand(lhs, 255_int32), amount), 255_int32)
  end function shl8

  integer(int32) function shr8(lhs, amount) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, amount
    value = ishft(iand(lhs, 255_int32), -amount)
  end function shr8

end program main
